#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright (c) 2011, 2012 Nick Semenkovich <semenko@alum.mit.edu> / WUSTL
# 
# Developed for the Gordon Lab, Washington University in St. Louis (WUSTL)
# http://gordonlab.wustl.edu/
#
# This software is released under the MIT License <http://www.opensource.org/licenses/mit-license.php>
#
"""
Some envelope stuff with mass spec stuff.
"""
from __future__ import absolute_import, division, print_function, with_statement # unicode_literals

__author__ = 'Nick Semenkovich <semenko@alum.mit.edu> and Gabriel Simon <gabrielmsimon@gmail.com>'
__copyright__ = 'Gordon Lab at Washington University in St. Louis / gordonlab.wustl.edu'
__license__ = 'MIT'
__version__ = '1.3.1'

from optparse import OptionParser, OptionGroup
from pkg_resources import parse_version
from xml.etree import ElementTree
import array
import base64
import csv
import datetime
import gc
import socket
import hashlib
import logging
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import multiprocessing
import os
import subprocess
import random
import re
import time
import signal
import string
import struct
import sys

### ---------------------------------------------
### Configuration Variables.
### You may need to edit this block for your setup.
### ---------------------------------------------

# http://physics.nist.gov/cuu/Constants/
MASS_PROTON = 1.007276466812

# http://www.nist.gov/pml/data/comp.cfm
# 14N = 14.0030740048
# 15N = 15.0001088982
N15_MASS_SHIFT = 0.997034893

AA_TO_N = {
    'A': 1,
    'R': 4,
    'N': 2,
    'D': 1,
    'C': 1,
    'E': 1,
    'Q': 2,
    'G': 1,
    'H': 3,
    'I': 1,
    'L': 1,
    'K': 2,
    'M': 1,
    'F': 1,
    'P': 1,
    'S': 1,
    'T': 1,
    'W': 2,
    'Y': 1,
    'V': 1,
    'U': 1,
    'O': 3,
}

# In daltons, a windows to choose MS1 scans
# e.g. a window of 20 on 500 da will be from [490:510]
MS1_WINDOW = 20


## *****
## Config Options
## *****

# isodist prediction range: These are the guesses we tell isodist to make regarding enrichment.
# These MUST have corresponding res_15Nshift_XXX.txt files in the /isodist/ folder!
#
# Node: You can add or remove values here, but execution (especially graph generation & parsing
# of fit files) will slow down substantially.
N_PERCENT_RANGE = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

# In Daltons, for guessing MS1 peak range from MS2 m/z
PEAK_SPREAD = 10

# Should we submit jobs to DRMAA? If false, we just run locally. This depends on DRMAA.
USE_DRMAA = False

# If a program has thread options (e.g. bowtie2 with --threads <int>), how many
# threads should we spawn?
# Of note -- Older versions of SGE may have issues with allocating >1 core/job.
MAX_THREADS = 1

### ---------------------------------------------
### Exceptions & Handlers
### ---------------------------------------------

def handle_SIGINT(signal, frame):
    """
    Caught a SIGINT / user pressed Ctrl-C
    """
    raise FatalError('User abort: Caught SIGINT')

class FatalError(Exception):
    """ Thrown when we should die. """
    global DRMAA_RUNNING_JOBS

    def __init__(self, message):
        Exception.__init__(self)

        fatal_logger = logging.getLogger('FatalError')

        if USE_DRMAA:
            # Terminate any running jobs.
            for job_id in DRMAA_RUNNING_JOBS:
                #noinspection PyUnresolvedReferences
                DRMAA_SESSION.control(job_id, drmaa.JobControlAction.TERMINATE)
                fatal_logger.critical('Killed SGE/DRMAA job %s' % job_id)

        fatal_logger.error(message)
        print('Terminating ...')
        sys.exit(1)

class ColorFormatter(logging.Formatter):
    """
    Colorize the logger class
    Via StackOverflow: http://stackoverflow.com/questions/384076/how-can-i-make-the-python-logging-output-to-be-colored
    """
    FORMAT = ("$BOLD%(name)-25s$RESET: %(levelname)-20s: "
              "%(message)-100s "
              "($BOLD%(filename)s$RESET:%(lineno)d)")
    
    #noinspection PyTupleAssignmentBalance
    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
    
    RESET_SEQ = "\033[0m"
    COLOR_SEQ = "\033[1;%dm"
    BOLD_SEQ = "\033[1m"
    
    COLORS = {
        'WARNING': YELLOW,
        'INFO': GREEN,
        'DEBUG': BLUE,
        'CRITICAL': YELLOW,
        'ERROR': RED
        }
    
    def formatter_msg(self, msg, use_color = True):
        """ Return a warning-based colorized message. """
        if use_color:
            msg = msg.replace("$RESET", self.RESET_SEQ).replace("$BOLD", self.BOLD_SEQ)
        else:
            msg = msg.replace("$RESET", "").replace("$BOLD", "")
        return msg

    def __init__(self, use_color=True):
        msg = self.formatter_msg(self.FORMAT, use_color)
        logging.Formatter.__init__(self, msg)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname
        if self.use_color and levelname in self.COLORS:
            fore_color = 30 + self.COLORS[levelname]
            levelname_color = self.COLOR_SEQ % fore_color + levelname + self.RESET_SEQ
            record.levelname = levelname_color
        return logging.Formatter.format(self, record)


### ---------------------------------------------
### Main program code.
### ---------------------------------------------

def main():
    """The main part of our MS analysis pipeline."""
    global USE_DRMAA, DRMAA_SESSION

    # Begin timing execution

    starttime = time.time()
    parser = OptionParser(usage="usage: %prog [options] input_directory\n\nInput directory must contain:\n\tDTASelect-filter.txt\n\t*.sqt\n\t*.mzXML",
                          version="%prog 1.0")

    parser.add_option("-v", "--verbose", help="Be verbose.",
                      default=logging.INFO, action="store_const", const=logging.DEBUG, dest="loglevel")

    # Enabling GC may save memory, but slows down computation.
    parser.add_option("--enable-gc", help="Enable garbage collection.",
                      default=False, action="store_true", dest="enable_gc")

    # This explicitly ignores hyperthreading pseudo-cores since isodist and matplotlib presumably hose the ALU.
    parser.add_option("--num-threads", help="Maximum number of children to spawn for parallelization. [Default: # of CPU cores]",
                      default=os.sysconf('SC_NPROCESSORS_ONLN'), action="store", type="int", dest="num_threads")

    group = OptionGroup(parser, "Skip Sections (optional)")
    group.add_option("--skip-isodist", help="Do not run isodist. (Isodist results must already exist.)",
                      default=False, action="store_true", dest="skip_isodist")
    group.add_option("--skip-graphs", help="Do not generate graphs.",
                      default=False, action="store_true", dest="skip_graphs")
    parser.add_option_group(group)

    # Parse the input and check.
    #noinspection PyTupleAssignmentBalance
    (options, args) = parser.parse_args()

    if not len(args) == 1:
        parser.error("You must specify exactly one input directory.\n\nTry --help for help.")

    input_directory = args[0]
    if not os.access(input_directory, os.R_OK):
        parser.error("Cannot read input directory.")

    # Check the #cpus requested
    if not 0 <= options.num_threads <= 128:
        parser.error("Invalid range for CPU limit. Choose from 0-128.")

    # Make sure we have some appropriately named input files, so we don't die later.
    directory_list = os.listdir(input_directory)
    if 'DTASelect-filter.txt' not in directory_list:
        parser.error("DTASelect-filter.txt not found in input directory.")
    if not len([fname for fname in directory_list if fname.endswith('.mzXML')]) == 1:
        parser.error("Exactly one .mzXML file must be present in the input directory.")
    if not len([fname for fname in directory_list if fname.endswith('.sqt')]) >= 1:
        parser.error(".sqt file(s) not found in input directory.")


    # Let's set up a logging system
    # We log DEBUG and higher to log file, and write INFO and higher to console.
    datestamp = datetime.datetime.now().strftime("%m%d-%H%M")
    logfilename = datestamp + '.' + socket.gethostname().split('.')[0]
    logging.basicConfig(filename = logfilename + '.log', filemode = 'w',
                        format = '%(asctime)s: %(name)-25s: %(levelname)-8s: %(message)s',
                        level = logging.DEBUG)

    # Define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    # level DEBUG if the user added --verbose or -v
    console.setLevel(options.loglevel)
    # Point to our colorized formatter.
    console.setFormatter(ColorFormatter())
    logging.getLogger('').addHandler(console)

    log_main = logging.getLogger('main')
    
    log_main.info('Welcome to enveloper.py!')
    log_main.info('Written by %s' % (__author__,))
    log_main.info('Developed for the %s' % (__copyright__,))
    log_main.info('Logging to %s' % (logfilename + '.log',))

    # I hate mid-program imports like this...
    if USE_DRMAA:
        import drmaa
        DRMAA_SESSION = drmaa.Session()
        DRMAA_SESSION.initialize()

    # Try to catch SIGINT for a more graceful shutdown
    signal.signal(signal.SIGINT, handle_SIGINT)

    # Check the version numbers, etc.
    pre_run_version_checks()

    
    ###### More application-specific functions
    # Parse DTA Select results file
    dta_select_data = parse_DTASelect(input_directory.rstrip('/') + "/DTASelect-filter.txt")

    if options.enable_gc:
        log_main.warning('Garbage collection is enabled, since --enable-gc was passed.')
        if sys.hexversion < 0x02070300:
            log_main.warning('Garbage collection is EXTREMELY slow on Python < 2.7')
    else:
        log_main.info('Disabling garbage is disabled. Use --enable-gc to override.')
        gc.disable()

    # Now that we have the isodist predicted spectra, parse the mzXML
    # TODO: Modularize these filenames better.
    ms1_data, ms2_to_ms1 = parse_mzXML(input_directory.rstrip('/') + '/' + [fname for fname in directory_list if fname.endswith('.mzXML')][0])

    # Now that we have the ms1 & DTASelect data, let's try to pick some peaks.
    # This is tricky, since DTASelect data is from MS2, so we have to kinda' guess the MS1 spectra.
    peptide_dict = extract_MS1_peaks(dta_select_data, ms1_data, ms2_to_ms1)

    log_main.debug('Enabling GC and marking ms1_data for deletion.')
    gc.enable() # Enable GC (if it was disabled), since the big stuff is done.
    del ms1_data # This is big. Go away. (Note: This doesn't imply python will free() the memory.)

    # Run isodist unless we're told to skip it. This can be slow.
    if options.skip_isodist:
        log_main.warning('Skipping isodist as requested. Assuming results already exist (may be stale).')
    else:
        run_isodist(dta_select_data, peptide_dict, options.num_threads)


    # Let's read those isodist results! WHOOO
    itime = time.time()
    isodist_results = read_isodist_results('./isodist/', peptide_dict, options.skip_graphs)
    log_main.info("reading isodist results took: %0.2f secs." % (time.time() - itime))

    # Why not make some matplotlib graphs?
    if options.skip_graphs:
        log_main.warning('Skipping peak graph generation as requested. (Graphs may be stale.)')
    else:
        make_peak_graphs(peptide_dict, isodist_results, options.num_threads)

    # Choose winners: rank predictions and choose the best FRC_NX value
    enrichment_predictions = pick_FRC_NX(peptide_dict, isodist_results)

    del isodist_results # This is huge.

    # Save output as CSV & HTML.
    generate_output(peptide_dict, enrichment_predictions, logfilename)

    # Cleanup.
    if USE_DRMAA:
        DRMAA_SESSION.exit()

    log_main.info("Execution took: %0.2f secs." % (time.time() - starttime))


def pre_run_version_checks():
    """
    Check that our dependencies are available and meet minimum version requirements.
    """
    log_prerun = logging.getLogger('pre_run_version_checks')

    # You should be able to use Python < 2.7.3, so long as you *disable* garbage collection.
    # In Python < 2.7, GC *kills* performance, probably 2/2 http://bugs.python.org/issue4074
    if sys.hexversion < 0x02070300:
        raise FatalError('Outdated Python version. Please use >=2.7.3')

    ### Tool Checks
    # Make sure isodist exists and runs.

    # Check isodist version (only one 2008 version at time of writing)
    env = dict(os.environ)
    env['LD_LIBRARY_PATH'] = './bin/'

    try:
        isodist_version_process = subprocess.Popen(['./bin/isodist_x86linux64'],
                                                   stdout=subprocess.PIPE,
                                                   stderr=subprocess.PIPE,
                                                   cwd='./isodist/',
                                                   env=env)
    except OSError:
        raise FatalError('isodist not found or not executable.')

    isodist_version_process.wait()
    stdout = isodist_version_process.communicate()[0]
    local_isodist_version = stdout[23:27] # Snip the isodist version string.

    if not cmp(parse_version(local_isodist_version), parse_version('2008')) >= 0:
        raise FatalError('isodist is outdated. Please use a version >= %s' % TOOLS_AND_VERSIONS['isodist'][1])

    # Make an isodist intermediate file directories
    for dirname in ['peaks', 'batch', 'input']:
        try:
            os.mkdir('./isodist/%s' % (dirname,))
        except OSError:
            if not os.access("./isodist/%s" % (dirname,), os.W_OK):
                raise FatalError('Unable to write to ./isodist/%s -- Make sure it exists and is writeable.' % (dirname,))

    # We shipped with two config files: exp_atom_defs.txt and res_15Nshift_XXX.txt.
    # The third required file (15Nshift_XXX.in) is dynamically generated.
    # Make sure they're unmodified, or warn users they've changed.
    isodist_files = [('./isodist/exp_atom_defs.txt', 'bbd69fd559741d93f0856ad6b9d7f8e8'),
                    ('./isodist/res_15Nshift_0.txt', 'c122efa100c61910dcfa7452415576c3'),
                    ('./isodist/res_15Nshift_20.txt', '6eabe6529ae1c972b6828065d34e3c99'),
                    ('./isodist/res_15Nshift_50.txt', '14e4ea1dac481dc4db1ba0a603376d74'),
                    ('./isodist/res_15Nshift_80.txt', 'c139deac216d13b6bf90f0041837fe1b'),
                    ('./isodist/res_15Nshift_100.txt', '67d4750db22afac837208bbc2c5a7da7'), ]
    isodist_observed_hashes = [hashlib.md5(file(fname).read()).hexdigest() for fname, _ in isodist_files]

    for input_file_pair, observed_hash in zip(isodist_files, isodist_observed_hashes):
        fname, valid_hash = input_file_pair
        if not valid_hash == observed_hash:
            # Are you seeing this error?
            # That's OK if you modified the input files. Otherwise, the are corrupt.
            log_prerun.warn('isodist input file has been modified: %s' % (fname,))


    log_prerun.debug('Version checks passed.')
    return True


def parse_DTASelect(DTASelect_file):
    """
    Open and parse a filtered DTA Select file.
      Input: path to a DTA elect-filtered file
      Output: A dict of the input file.
    """
    parse_dta_log = logging.getLogger('parse_DTASelect')
    parse_dta_log.info('Parsing DTA Select file: %s' % (DTASelect_file,))

    # Store a dictionary of the DTASelect-filter file
    # This is structured as:
    # protein_id: {metadata: , peptide:}
    #      meta:
    dta_dict = {}

    dta_select_csv = csv.reader(open(DTASelect_file, 'rb'), delimiter = '\t')

    # Temporary parsing variables
    past_header = False
    current_keys = []
    added_peptides = False
    peptide_dict = {}

    # These files aren't the easiest to parse. We pick lines based on the number of TSV elements.
    for line in dta_select_csv:
        # Have we gone past the header in the file yet?
        if past_header:
            if len(line) == 9:
                # Length 9 lines are protein groups (they delimit peptide sections)
                if not added_peptides:
                    # This section has multiple protein headers, or we've just started to parse the file.
                    current_keys.append(line[0])
                else:
                    # We must've just entered a new protein section. Reset our "current_keys" list.
                    added_peptides = False
                    for key in current_keys:
                        dta_dict[key]['peptides'] = peptide_dict
                    current_keys = [line[0]]
                    peptide_dict = {}

                if line[0] in dta_dict:
                    # I don't think this is ever possible. But let's be paranoid.
                    raise FatalError('Duplicate protein key in DTA Select file!')

                # Make a dict in our dict!
                dta_dict[line[0]] = {'metadata': {}, 'peptides': {}}

                # Set our metadata. In the file, the format is:
                # 'Sequence Count', 'Spectrum Count', 'Sequence Coverage', 'Length', 'MolWt', 'pI', 'Validation Status', 'Descriptive Name'
                metadata_dict = dta_dict[line[0]]['metadata']
                keys = ['seq_count', 'spect_count', 'seq_cov', 'length', 'molwt', 'pI', 'validated', 'name']
                types = [int, int, str, int, int, float, str, str]
                for key, value, cast in zip(keys, line[1:], types):
                    metadata_dict[key] = cast(value)
 

            elif len(line) == 13:
                # Length 13 lines are peptide data, which we add to a protein entry.

                # Mark that we've started adding data. Otherwise we'd add it to the wrong key!
                added_peptides = True

                # In the file, the format is:
                # The FileName is our unique dict key.
                # ['Unique', 'FileName', 'XCorr', 'DeltCN', 'Conf%', 'M+H+', 'CalcM+H+', 'TotalIntensity', 'SpR', 'Prob Score', 'IonProportion', 'Redundancy', 'Sequence']
                # TODO: Restructure this for performance. This is extremely un-Pythonic.
                peptide_key = line[1]
                del line[1] # This is not friendly.

                # Don't think this is possible, but let's be paranoid.
                if peptide_key in peptide_dict:
                    raise FatalError('Duplicate FileName key in DTASelect-filter.txt')

                keys = ['unique', 'xcorr', 'delt_cn', 'conf', 'mh', 'calc_mh', 'tot_intensity', 'spr', 'prob_score', 'ion_proportion', 'redundancy', 'sequence']
                types = [bool, float, float, float, float, float, float, float, float, float, int, str]
                peptide_dict[peptide_key] = {}
                for key, value, cast in zip(keys, line, types):
                    peptide_dict[peptide_key][key] = cast(value)

            elif len(line) == 4:
                # We're at the end of the file. Victory is ours!
                for key in current_keys:
                    dta_dict[key]['peptides'] = peptide_dict
                break
            else:
                raise FatalError('Odd structure in DTA Select file!')
                
        # We aren't entirely though the header yet.
        else:
            # The header is 13 elements and starts with "Unique"
            if len(line) == 13 and line[0] == "Unique":
                past_header = True
                parse_dta_log.debug('Found data stream in DTA Select file')


    parse_dta_log.info('Finished parsing')

    return dta_dict

def extract_MS1_peaks(dta_select_data, ms1_data, ms2_to_ms1):
    """
    Given the dta_select data, and the ms1 from the mzXML,
    try to extract representative peaks for each peptide.
    """
    extract_peak_logger = logging.getLogger('extract_MS1_peaks')
    peptide_dict = {}
    
    az_only_pattern = re.compile('[^A-Z]+')

    # Calculate some charge distributions
    charge_dist = [0]*6

    # We don't need the protein key, yet. Maybe in a future version.
    # pylint: disable=W0612
    for protein_key, protein_data in dta_select_data.iteritems():
        for peptide_key, peptide_data in protein_data['peptides'].iteritems():

            # Skip already done peptides:
            # These exist when there are >1 predicted *proteins* in the DTASelect-file
            if peptide_key in peptide_dict:
                extract_peak_logger.debug('Peptide from %s already predicted. Skipping.' % (peptide_key,))
                continue

            # We split the key and discard the .sqt filename
            #noinspection PyTupleAssignmentBalance
            scan_start, scan_stop, charge = [int(elt) for elt in peptide_key.split('.')[1:]]
            if scan_start != scan_stop:
                raise FatalError('Unsupported MS2 scan range found in DTASelect')
            # Since we don't always have a 7:1 ms2:ms1 ratio, we can't do this awesome modulo. :(
            # parent_scan = scan_start - ((scan_start - 1) % 8)
            parent_scan = ms2_to_ms1[scan_start]
            extract_peak_logger.debug('Parent: %s for MS2: %s' % (parent_scan, scan_start))

            if 1 <= charge <= 5:
                charge_dist[charge] += 1
            else:
                extract_peak_logger.warning('Charge of %s seen in MS2 scan %s' % (charge, scan_start))

            calc_mh = peptide_data['calc_mh'] # This is the +1 state from the DTASelect file
            # Strip protease cleavage sites and non-[A-Z] characters from the sequence
            peptide_sequence = az_only_pattern.sub('', peptide_data['sequence'][2:-2])
            calc_mz = (calc_mh + ((charge - 1) * MASS_PROTON)) / charge
            n_in_peptide = sum([AA_TO_N[aa] for aa in peptide_sequence])
            n15_adjusted_mass = calc_mz + (n_in_peptide * N15_MASS_SHIFT)

            # Print for debugging
            extract_peak_logger.debug('\t CalcMH: %s (Seq: %s, Charge: %s)' % (calc_mh, peptide_sequence, charge))
            extract_peak_logger.debug('\t M/Z: %s' % (calc_mz,))
            extract_peak_logger.debug('\t 15N M/Z: %s' % (n15_adjusted_mass,))

            # Now let's pick a range of probable MS1 peaks:
            try:
                mz_from_parent = ms1_data[parent_scan]['peak']
            except KeyError:
                raise FatalError('Parent scan not in MS1 dict: Parent calculation wrong?')

            # Take from w/i our range
            extracted_ms1 = [(m, z) for m, z
                                    in mz_from_parent
                                    if (calc_mz - MS1_WINDOW/2.0) < m < (n15_adjusted_mass + MS1_WINDOW/2.0)]

            peptide_dict[peptide_key] = {'sequence': peptide_sequence,
                                         'mz': calc_mz,
                                         'n15mz': n15_adjusted_mass,
                                         'charge': charge,
                                         'peaks': extracted_ms1, }

    extract_peak_logger.info('Summary of charge distribution:')
    for index, item in enumerate(charge_dist):
        extract_peak_logger.info('\t%s: %0.2f%%' % (index, (item / float(sum(charge_dist)) * 100)))

    return peptide_dict

def make_peak_graphs(peptide_dict, isodist_results, num_threads):
    """
    Make some graphs of the peaks. Now featuring parallelism!

    Note: These threads take a *lot* of RAM. Like 8gb/thread.

    FYI: Threadsafety here is not 100% certain. This may break at high thread numbers (>=24).
      Although I think that bug was fixed in Python >=2.7.
    """

    graph_log = logging.getLogger('make_peak_graphs')
    graph_log.info('Generating graphs. This may take a few minutes.')

    graph_log.debug('Making output directories in /graphs/.')
    # Let's create output directories for graphs
    for peptide_key in peptide_dict.iterkeys():
        try:
            # Makedirs will make the parent /graphs/ if it doesn't exist.
            os.makedirs('./graphs/%s' % (peptide_key))
        except OSError:
            # Dir might already exist. If it's unwriteable, we'll FATAL it later.
            pass

    pool = multiprocessing.Pool(num_threads)
    tasks = [(key, val, isodist_results[key]) for key, val in peptide_dict.iteritems()]
    results = []

    chunk_size = len(tasks)//num_threads
    graph_log.debug('Task size is %s with %s threads' % (len(tasks), num_threads))
    graph_log.debug('Using chunk_size of %s' % (chunk_size,))

    # TODO: Use error_callback?
    r = pool.map_async(_peak_graph_cmd, tasks, chunk_size, callback=results.append)
    r.wait() # Block until our pool returns

    try:
        if len(results[0]) != len(tasks):
            # You could take a set intersection and see what didn't return.
            raise FatalError('A graphing thread failed!')
    except IndexError:
        raise FatalError('A graphing thread failed: Make sure the ./graphs/ output directories are writeable.')
 
    graph_log.info('Graphs generated successfully.')
    
    return True

def pick_FRC_NX(peptide_dict, isodist_results):
    """
    Choose the best isodist FRC_NX enrichment guess. There isn't one clean way to do this.
    I opt for:
     Rank by chi_sq from isodist
     Ensure top four predictions are within 1%
       If so, choose the average of the four.
       If not, refuse to make an FRC_NX guess.
    """

    frc_nx_log = logging.getLogger('pick_FRC_NX')
    frc_nx_log.info('Determining optimal FRC_NX enrichment percentages.')

    tasks = [(key, val, isodist_results[key]) for key, val in peptide_dict.iteritems()]
    
    frc_nx_log.info('Predictions were made for %s peptides.' % (len(tasks),))

    enrichment_predictions = {}

    # Loop over the tasks based on key (key: g08_L50.21917.21917.2)
    #
    # I'm not sure what the best way to settle on a prediction is. The options include:
    #  - Direct comparsions to raw peaks (ignoring the isodist CHI_SQ scores)
    #  - Averaging/other statistics based on CHI_SQ or other values
    #  - Windowing/mode selection
    #
    # The first option never worked great for me, so I've settled on using a windowing-type approach
    # followed by an average. We perform the following algorithm:
    #  1. For the 10 N_PERCENT_RANGE values [0, 10, 20 ...], take each predicted FRC_NX and ...
    #  2. Using the current FRC_NX, see if the next FRC_NX is within 1%:
    #    - If yes, continue until >= 5 values are within 1% of each other
    #    - If no, move to the next FRC_NX and go back to step #2
    #  3. If >=5 FRC_NX predictions are within 1% of each other, return their average, otherwise,
    #    we refuse to maek an FRC_NX enrichment prediction.
    #
    # This may seem inelegant, but in practice, is is extremely stringent.
    for k, _, i in tasks:
        frc_nx_log.debug('Choosing enrichment for %s' % (k,))
        # Set key of peptide id
        enrichment_predictions[k] = {}
        # Add raw FRC_NX guesses from isodist
        for percent in N_PERCENT_RANGE:
            enrichment_predictions[k]['guess_' + str(percent)] = i[percent]['frc_nx']

        # Keep some stats for std dev
        sum_all_enrich = 0.0
        sum_all_enrich_sq = 0.0
        
        sum_window_sq = 0.0

        # Loop over the N_PERCENT_RANGE isodist guesses with our sliding window 1% limit
        current_window_percent = -10  # An impossible to return enrichment value.
        enrichment_window = []
        for percent in N_PERCENT_RANGE:
            sum_all_enrich += i[percent]['frc_nx']
            sum_all_enrich_sq += i[percent]['frc_nx'] ** 2

            if (current_window_percent - .01 ) < i[percent]['frc_nx'] < (current_window_percent + .01):
                # We're within 1% of the previous value
                enrichment_window.append(i[percent]['frc_nx'])
                current_window_percent = i[percent]['frc_nx']
                sum_window_sq += i[percent]['frc_nx'] ** 2
            elif (len(enrichment_window) >= 5):
                # We have five or more values that are awesome. We're set. Let's get outa' here.
                break
            else:
                # We're not within 1%. Clear the window. Let's hope there are more values!
                enrichment_window = []
                current_window_percent = i[percent]['frc_nx']
                sum_window_sq = 0.0

        # Ok, let's see what that loop left in our enrichment window ...
        frc_nx_log.debug('\tWindow: %s items, %s ' % (len(enrichment_window), enrichment_window))
        frc_nx_log.debug('\tRaw: %s' % ([i[percent]['frc_nx'] for percent in N_PERCENT_RANGE]))

        if len(enrichment_window) >= 5:
            mean = sum(enrichment_window)/float(len(enrichment_window))
            enrichment_predictions[k]['percent'] = mean
            frc_nx_log.debug('\tChoosing: %0.2f%%' % (mean * 100,))
            frc_nx_log.debug('\t\tStd Dev: %0.5f' % (math.sqrt((sum_window_sq / len(enrichment_window)) - (mean ** 2)))) 
        else:
            all_mean = sum_all_enrich/float(len(N_PERCENT_RANGE))
            frc_nx_log.warn('\tPrediction failed for %s' % (k,))
            frc_nx_log.debug('\t\tOverall mean: %0.2f' % (all_mean,))
            frc_nx_log.debug('\t\tStd Dev: %0.2f' % (math.sqrt((sum_all_enrich_sq / len(N_PERCENT_RANGE)) - (all_mean ** 2))))

    frc_nx_log.info('Percentages chosen successfully.')
    
    return enrichment_predictions

def generate_output(peptide_dict, enrichment_predictions, logfilename):
    """
    Save a final output summary of our predictions. This outputs as CSV & HTML.
    """

    out_csv_fname = 'results/%s/%s' % (logfilename, 'envelope_results.csv')
    out_html_fname = 'results/%s/%s' % (logfilename, 'envelope_results.html')

    output_log = logging.getLogger('generate_output')
    output_log.info('Saving CSV output as: results/%s' % (out_csv_fname,))

    # Make a results directory
    try:
        os.mkdir('./results/%s' % (logfilename,))
    except OSError:
        pass


    # Prepare a dict for writing.
    # Add a 'key' value to the dict for our DictWriter
    for key in peptide_dict.iterkeys():
        peptide_dict[key]['id'] = key

        # Merge in any of the enrichment_predictions, overwrite collisions (!)
        for enrich_key, enrich_val in enrichment_predictions[key].iteritems():
            peptide_dict[key][enrich_key] = enrich_val


    # Dictionary keys for our output results.
    output_keys = ['id', 'sequence', 'charge', 'mz', 'n15mz', 'percent']
    # Add the "guess_0" -> 100 columns
    for percent in N_PERCENT_RANGE:
        output_keys.append('guess_' + str(percent))

    # Open & write our CSV file
    # Protip: Write to sys.stderr if you're debugging this.
    with open(out_csv_fname, 'wb') as csvout:
        csv_out = csv.DictWriter(csvout, output_keys, extrasaction='ignore')
        # Here's a header
        csv_out.writeheader()
        # And the contents
        csv_out.writerows(peptide_dict.itervalues())
    output_log.info('CSV sucessfully generated.')


    # Let's make some HTML, too.
    # Note: This uses DataTables -- if the number of columns changes, the table may break.
    #   Be sure to update the HTML headers if you change things.
    output_log.info('Saving HTML output as: %s' % (out_html_fname,))
    with open(out_html_fname, 'wb') as htmlout, open('.html/_html_header.html', 'rb') as headers, open('.html/_html_footer.html', 'rb') as footers:
        htmlout.writelines(headers.readlines())
        for elt in peptide_dict.itervalues():
            print('<tr><td>', file=htmlout)
            print('</td><td>'.join([str(elt.get(x, '')) for x in output_keys]), file=htmlout)
            print('</td></tr>', file=htmlout)
        htmlout.writelines(footers.readlines())
    output_log.info('HTML successfully generated.')


    return True

def _peak_graph_cmd(task):
    """
    Generate peak graphs as part of the multiprocessing/map_async command
    """
    peptide_key, peptide_value, isodist_results = task

    graph_thread_log = logging.getLogger('_peak_graphs_cmd')
    graph_thread_log.debug('Graphing %s' % peptide_key)
    
    # Not that this app is secure, but just so we don't break things, sanitize some outputs.
    valid_chars = "-_.()%s%s" % (string.ascii_letters, string.digits)

    for n_percent in N_PERCENT_RANGE:
        # Set up our figure
        fig = matplotlib.pyplot.figure(figsize=(12, 5))
        fig.suptitle('%s' % (peptide_key, ))
        matplotlib.pyplot.xlabel('M/Z')
        matplotlib.pyplot.ylabel('Intens')

        # Add some metadata in a box
        matplotlib.pyplot.annotate("Seq: %s\nMZ: %0.4f\n15N MZ: %0.4f" %
                                   (peptide_value['sequence'], peptide_value['mz'], peptide_value['n15mz']),
            (0.85, 0.95),
            xycoords="axes fraction", va="center", ha="left",
            bbox=dict(boxstyle="round, pad=1", fc="w"))

        # Let's also grab the metadata from isodist
        isodist_data = isodist_results[n_percent]

        matplotlib.pyplot.annotate("FRC_NX: %0.4f\nCHI Sq: %e\nAMP_U: %0.4f\nAMP_L: %0.4f" %
                                   (isodist_data['frc_nx'], isodist_data['chisq'], isodist_data['amp_u'], isodist_data['amp_l']),
            (0.85, 0.70),
            xycoords="axes fraction", va="center", ha="left",
            bbox=dict(boxstyle="round, pad=1", fc="w"))

        # Add our actual data (we need to add a Plot to our Figure)
        ax = fig.add_subplot(111)
        #noinspection PyTupleAssignmentBalance
        m, z = zip(*peptide_value['peaks'])
        ax.plot(m, z, 'r-', linewidth = 0.5)

        # Get the m/z data.
        # Note: These are pre-filtered to exclude all pairs with intensity <=0, as they mess w/ the graph.
        isodist_m = [m/float(peptide_value['charge']) for m in isodist_data['peak_fit'][0]]
        isodist_z = isodist_data['peak_fit'][1]

        ax.plot(isodist_m, isodist_z, 'b-.', linewidth = 1.2)
        #ax.autoscale_view() # I really don't know what this does.
        ax.grid(True)

        # Semi-sanitized. Not secure, but we aren't likely to accidentally break things if peptide_key is insane.
        safer_key = ''.join(c for c in peptide_key if c in valid_chars)

        matplotlib.pyplot.savefig("graphs/%s/%s.png" % (safer_key, n_percent))

        # I'm not sure if this matters for memory usage.
        # Matplotlib uses a *lot* of RAM and nothing seems to reduce it.
        matplotlib.pyplot.close('all')

    return True

def run_isodist(dta_select_data, peptide_dict, num_threads):
    """
    Run isodist.
    """
    isodist_log = logging.getLogger('run_isodist')
    isodist_log.info('Running isodist')

    # Find the unique peptides from the dta_select_data dict
    # TODO: Make this a little easier to follow when accessing data?
    unique_peptides = set()
    for data in dta_select_data.values():
        for peptide_data in data['peptides'].values():
            unique_peptides.add(peptide_data['sequence'][2:-2])

    isodist_log.info('%s unique peptide sequences found' % (len(unique_peptides),))

    # Isodist settings
    isodist_settings = {'number_of_iterations': 10,
                        'guess_for_baseline': "150.0 auto",
                        'accuracy_offset': "0.01 variable",
                        'gaussian_width': "0.003 variable",
                        'n_percent': 0000, } # Note: n_percent is adjusted on-the-fly below.

    input_template = "fitit\nbatch/%(batchfile_name)s/%(n_percent)s.batch\nexp_atom_defs.txt\n" + \
                     "res_15Nshift_%(n_percent)s.txt\n%(number_of_iterations)s\n100.0\n%(guess_for_baseline)s\n" + \
                     "%(accuracy_offset)s\n%(gaussian_width)s"

    isodist_log.info('Creating isodist input files in ./isodist/')
    for peptide_key, peptide_value in peptide_dict.iteritems():
        for dirname in ['peaks', 'batch', 'input']:
            try:
                os.mkdir('./isodist/%s/%s' % (dirname, peptide_key))
            except OSError:
                # Dir might already exist. If it's unwriteable, we'll FATAL it below.
                pass

        # We need multiple copies of these files -- one pair per estimated N%.
        # TODO: Sanitize these values.
        try:
            for n_percent in N_PERCENT_RANGE:
                peaks = open('./isodist/peaks/%s/%s.tsv' % (peptide_key, n_percent), 'w')
                batch = open('./isodist/batch/%s/%s.batch' % (peptide_key, n_percent), 'w')
                isoinput = open('./isodist/input/%s/%s.in' % (peptide_key, n_percent), 'w')

                # We could rewrite isodist to take more flags, but for now let's just use the default (messy) format.
                # We make so many input files so we have discrete work units for parallel computation.

                # Write the raw peaks -- we only need one of these files.
                for m, z in peptide_value['peaks']:
                    print("%s\t%s" % (m, z), file=peaks)

                # Make the .in file, which is a big, annoying template
                isodist_settings['batchfile_name'] = peptide_key
                isodist_settings['n_percent'] = n_percent
                print(input_template % isodist_settings, file=isoinput)
                
                # Make the batchfile
                print("%s %s peaks/%s/%s.tsv" % (peptide_value['sequence'], peptide_value['charge'], peptide_key, n_percent), file=batch)

                peaks.close()
                batch.close()
                isoinput.close()
        except IOError:
            raise FatalError('Could not write an ./isodist/ file!')


    # Should we spawn jobs via DRMAA?
    if USE_DRMAA:
        isodist_log.info('Distributing isodist jobs via DRMAA/SGE.')
        raise FatalError('Not implemented.')
    else:
        isodist_log.info('Running isodist jobs locally, as DRMAA/SGE is disabled.')
        
        pool = multiprocessing.Pool(num_threads)
        #tasks = peptide_dict.keys()
        tasks = [x + "/" + str(y) for x in peptide_dict.keys() for y in N_PERCENT_RANGE]
        results = []

        # TODO: Use error_callback?)
        r = pool.map_async(_isodist_cmd, tasks, num_threads, callback=results.append)
        r.wait() # Block until our pool returns
        if len(results[0]) != len(tasks):
            # You could take a set intersection and see what didn't return.
            raise FatalError('isodist execution failed!')

        isodist_log.info('isodist executions complete')

    return True


def read_isodist_results(input_path, peptide_dict, skip_graphs):
    """
    Open our raw isodist results files and parse their delightful results.

    I've tried to speed this up as much as possible, but the peak .fit files simply take a while to read.
    """
    read_logger = logging.getLogger('read_isodist_results')
    read_logger.info('Reading isodist results. This will take a few minutes.')

    # Save our results in a dict
    isodist_results = {}

    # We'll save some of the .batch.csv results
    keys = ['mw', 'z', 'tim', 'chisq', 'symb', 'mz', 'b', 'off', 'gw', 'amp_u', 'amp_l', 'frc_nx']
    casts = [float, int, str, float, int, float, float, float, float, float, float, float]
    # We ran isodist on each peptide

    read_logger.debug('Reading %s percent values over %s keys.' % (len(N_PERCENT_RANGE), len(peptide_dict)))
    progress = 0
    total = len(N_PERCENT_RANGE)*len(peptide_dict)

    for peptide_key, n_percent in [(x, y) for x in peptide_dict.iterkeys() for y in N_PERCENT_RANGE]:
        progress +=1
        read_logger.debug('%0.2f%% complete. On key: %s, %s%%' % (progress/total, peptide_key, n_percent))

        # Make the [peptide_key] dict, if it doesn't exist.
        isodist_results.setdefault(peptide_key, {})

        isodist_results[peptide_key][n_percent] = {}

        # TODO: Change this to a directory tree?
        with open(input_path + 'batch/' + peptide_key + '/' + str(n_percent) + '.batch.csv', 'rb') as stats:
            statline = stats.readlines()[1]
            # This is parsing the .batch.csv metadata results.
            for key, value, cast in zip(keys, statline.split(',')[3:], casts):
                isodist_results[peptide_key][n_percent][key] = cast(value)

        # If we're skipping graphing, we don't need the enormous peak_fit data.
        # NOTE: This is VERY SLOW. Sorry.
        if not skip_graphs:
            # Now let's get the raw peak fit data, and add it to the dict, too!
            with open(input_path + 'peaks/' + peptide_key + '/' + str(n_percent) + '.fit', 'rb') as fit_file:
                peak_fit = csv.reader(fit_file, quoting=csv.QUOTE_NONNUMERIC) # Pre-cast to float
                # Note: We filter to intensity > 0, otherwise the graphs look strange.
                m, z = zip(*[(m, z) for m, z in peak_fit if z > 0])
                filtered_peaks_m = array.array('f', m)
                filtered_peaks_z = array.array('f', z)

            isodist_results[peptide_key][n_percent]['peak_fit'] = (filtered_peaks_m, filtered_peaks_z)


    read_logger.info('isodist results loaded successfully.')

    return isodist_results

        
def _isodist_cmd(infile):
    """
    Run isodist as part of the multiprocessing/map_async command
    NOTE: This is included here, as map_async/map handle nested functions poorly.
    """

    # Unspecified excepts are ugly, but we're paranoid here, because the parent thread isn't
    # the best about catching errors in isodist execution.
    isodist_cmd_log = logging.getLogger('_isodist_cmd')
    # Is this (& the instantiation) slow?

    # There's a library issue with isodist :(
    env = dict(os.environ)
    env['LD_LIBRARY_PATH'] = './bin/'

    isodist_cmd_log.debug("Running on %s" % (infile,))
    try:
        p = subprocess.Popen(['nice', './bin/isodist_x86linux64', './input/%s.in' % (infile,)],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             cwd='./isodist/',
                             env=env,
                             )
        stdoutdata, stderrdata = p.communicate() # This will block until isodist finishes.
    except:
        raise FatalError('isodist execution failed on %s' % (infile,))

    if p.returncode or stderrdata or ("FINAL FIT PARAMETERS" not in stdoutdata):
        raise FatalError('isodist failed on %s' % (infile,))

    return infile


def parse_mzXML(mzXML_file):
    """
    Open and parse an mzXML file.
    """
    parse_mzXML_log = logging.getLogger('parse_mzXML')
    parse_mzXML_log.info('Parsing mzXML file %s' % (mzXML_file,))

    # We'll return this dict of ms1 data at the end
    ms1 = {}
    # We also return an ms2->ms1 parent mapping
    ms2_to_ms1 = {}

    # Establish some dict key names & castings
    # Got rid of some "fixed" keys that never change, per:
    # http://sashimi.sourceforge.net/schema_revision/mzXML_3.2/mzXML_3.2.xsd
    ms1keys = ['polarity', 'basePeakIntensity', 'retentionTime', 'basePeakMz', 'peaksCount', 'lowMz', 'scanEvent', 'totIonCurrent', 'highMz', 'centroided']
    ms1types = [str, float, str, float, int, float, int, float, float, int]

    # Open the mzXML
    namespace = "{http://sashimi.sourceforge.net/schema_revision/mzXML_3.2}"
    parsed_xml = ElementTree.parse(mzXML_file)
    scans = parsed_xml.findall(".//%sscan" % namespace)

    parse_mzXML_log.info('Found %s scans in mzXML file' % (len(scans),))
    # Loop over each XML scan entry
    for scan in scans:
        scan_num = int(scan.attrib['num'])
        if not scan_num % 2000:
            parse_mzXML_log.debug('On scan %s' % (scan_num,))
        # Are we in a MS1 or MS2 file?
        if scan.attrib['msInstrumentID'] == "IC1":
            # Store our MS1 values in a dict for this scan
            ms1[scan_num] = {}
            # Cast & store our ms1 values
            for key, cast in zip(ms1keys, ms1types):
                ms1[scan_num][key] = cast(scan.attrib[key])

            # MS1 scans have only one "peak" child ("scan" is an iterable)
            # We check that these are uncompressed & 32-bit IEEE-754 network-order b64
            if int(scan[0].attrib['compressedLen']):
                raise FatalError('Sorry, compressed mzXML parsing is not implemented')
            if int(scan[0].attrib['precision']) != 32:
                raise FatalError('Sorry, 64-precision IEEE-754 parsing is not implemented.')
            if scan[0].attrib['byteOrder'] != "network":
                raise FatalError('Sorry, non-big-endian storage is not implemented.')
            if scan[0].attrib['pairOrder'] != "m/z-int":
                raise FatalError('Sorry, non m/z-int order is not implemented.')

            decoded_b64 = base64.b64decode(scan[0].text) # Decode the packed b64 raw peaks
            # These are packed as big-endian IEEE 754 binary32
            floating_tuple = struct.unpack('>' + str(len(decoded_b64)//4) + 'f', decoded_b64) 
            ms1[scan_num]['peak'] = zip(floating_tuple[::2], floating_tuple[1::2])

        elif scan.attrib['msInstrumentID'] == "IC2":
            # Store ms2 -> parent ms1 scan data. That's it.
            ms2_to_ms1[scan_num] = int(scan[0].attrib['precursorScanNum'])
        else:
            raise FatalError('Unknown msInstrumentID in mzXML.')

    parse_mzXML_log.info('Extracted %s ms1 scans.' % (len(ms1),))

    return ms1, ms2_to_ms1


def deploy_drmaa_job(job_command, job_parameters):
    """
    Submit these jobs a DRMAA commands. WARNING: Not currently implemented.
    If you needed expansive parallelism across machines, you could break out the isodist commands.
    """
    global DRMAA_RUNNING_JOBS

    drmaa_logger = logging.getLogger('deploy_drmaa_job')

    assert(os.access(job_command, os.X_OK))
    assert(type(job_parameters) is list)

    # We have a random delay so SGE has a chance get its emotions in check.
    # Not sure if this matters for the DRMAA interface, but hey, it's just a few seconds.
    time.sleep(random.random() / 2.0)
    
    job_template = DRMAA_SESSION.createJobTemplate()
    job_template.remoteCommand = job_command
    job_template.args = job_parameters
    
    job_id = DRMAA_SESSION.runJob(job_template)
    drmaa_logger.info("Deployed DRMAA job: %s" % job_id)
    DRMAA_RUNNING_JOBS.append(job_id)

    # Use runBulkJobs when parallelizing?
    DRMAA_SESSION.deleteJobTemplate(job_template)
    return job_id

### ---------------------------------------------

if __name__ == '__main__':
    # Let's establish some global variables to hold SGE/DRMAA runners and whatnots.
    CMD_LINE_INPUT_FILES = None
    if USE_DRMAA:
        # Store DRMAA job IDs in case we need to destroy them!
        DRMAA_RUNNING_JOBS = []
        # Global variable to hold a session. This is a little dirty.
        DRMAA_SESSION = None

    main()
