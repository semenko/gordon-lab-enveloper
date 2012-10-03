#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2011, 2012 Nick Semenkovich <semenko@alum.mit.edu> / WUSTL
#
# Developed for the Gordon Lab, Washington University in St. Louis (WUSTL)
# http://gordonlab.wustl.edu/
#
# This software is released under the MIT License:
#  <http://www.opensource.org/licenses/mit-license.php>
#
"""
Perform MS enrichment predictions on large datasets.
"""
from __future__ import absolute_import, division, print_function, with_statement

__author__ = 'Nick Semenkovich <semenko@alum.mit.edu> and Gabriel Simon <gabrielmsimon@gmail.com>'
__copyright__ = 'Gordon Lab at Washington University in St. Louis / gordonlab.wustl.edu'
__license__ = 'MIT'
__version__ = '1.3.1'

from base64 import b64decode
from optparse import OptionParser, OptionGroup
from struct import unpack
from xml.etree import cElementTree
import array
import csv
import cPickle
import datetime
import gc
import hashlib
import logging
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import multiprocessing
import os
import socket
import subprocess
import random
import re
import time
import signal
import shutil
import string
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

    BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

    RESET_SEQ = "\033[0m"
    COLOR_SEQ = "\033[1;%dm"
    BOLD_SEQ = "\033[1m"

    COLORS = {'WARNING': YELLOW,
              'INFO': GREEN,
              'DEBUG': BLUE,
              'CRITICAL': YELLOW,
              'ERROR': RED,
              }

    def formatter_msg(self, msg, use_color=True):
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

    # Some files are very slow to read (lots of tiny files), so we'd prefer a pickle cache.
    # This is useful if you find yourself usign the "--skip-isodist" flag
    #   (e.g. you're trying to modify the output HTML/CSV results.)
    parser.add_option("--force-cache", help="Cache and re-use parsed data. Useful if you are tweaking the HTML output.",
                      default=False, action="store_true", dest="force_cache")

    # Parse the input and check.
    #noinspection PyTupleAssignmentBalance
    (options, args) = parser.parse_args()

    if not len(args) == 1:
        parser.error("You must specify exactly one input directory.\n\nTry --help for help.")

    input_directory = args[0]
    if not os.access(input_directory, os.R_OK):
        parser.error("Cannot read input directory.")

    # Sanity check the #cpus requested
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
    results_path = datestamp + '.' + socket.gethostname().split('.')[0]
    logging.basicConfig(filename=results_path + '.log', filemode='w',
                        format='%(asctime)s: %(name)-25s: %(levelname)-8s: %(message)s',
                        level=logging.DEBUG)

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
    log_main.info('Logging to %s' % (results_path + '.log',))

    # Make a results directory early on, so we don't die at the last possible step.
    try:
        os.mkdir('./results/%s' % (results_path,))
    except OSError:
        parser.error('Results directory not writeable or already exists.')

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

    # We can short-circuit the mzXML parsing and peak extraction if we're super cached.
    peptide_dict = None
    if options.force_cache:
        # NOTE: Caching is essentially *only* useful if you're modifying this script.
        log_main.warning('You have requested to use cached data! Are you sure you want that?!')
        peptide_cache_key = hashlib.md5(input_directory).hexdigest() + '.peptide_dict.pkl'
        try:
            with open('.cache/' + peptide_cache_key, 'rb') as peptide_cached:
                peptide_dict = cPickle.load(peptide_cached)
        except IOError:
            log_main.info('No peptide_dict cache found. We\'ll write one momentarily.')

    # We almost always enter this chunk. The exception is if we've cached results,
    # probably for development purposes.
    if peptide_dict is None:
        # Parze the input mzXML data (we need this for peptide dict generation)
        # TODO: Modularize these filenames better.
        ms1_data, ms2_to_ms1 = parse_mzXML(input_directory.rstrip('/') + '/' + [fname for fname in directory_list if fname.endswith('.mzXML')][0])

        # Now that we have the ms1 & DTASelect data, let's try to pick some peaks.
        # This is tricky, since DTASelect data is from MS2, so we have to kinda' guess the MS1 spectra.
        peptide_dict = extract_MS1_peaks(dta_select_data, ms1_data, ms2_to_ms1)

        log_main.debug('Enabling GC and marking ms1_data for deletion.')
        gc.enable()   # Enable GC (if it was disabled), since the big stuff is done.
        del ms1_data  # This is big. Go away. (Note: This doesn't imply python will free() the memory.)

        if options.force_cache:
            with open('.cache/' + peptide_cache_key, 'wb') as pep_pickle_cache:
                cPickle.dump(peptide_dict, pep_pickle_cache)
    else:
        # Whoa, we got a peptide cache! Cool.
        log_main.warning('Using cached peptide data: If you\'ve modified input data, the cache will be stale!')

    # Run isodist unless we're told to skip it. This can be slow.
    if options.skip_isodist:
        log_main.warning('Skipping isodist as requested. Assuming results already exist (may be stale).')
    else:
        run_isodist(dta_select_data, peptide_dict, options.num_threads)

    # Let's read those isodist results! (Perhaps from a cache.)
    # Pickled' caches can be large, so we don't always write them.
    itime = time.time()
    isodist_results = False
    if options.force_cache:
        try:
            suffix = '.pkl'
            if options.skip_graphs:
                suffix = '.nographs.pkl'
            cache_key = hashlib.md5(input_directory).hexdigest() + suffix
            with open('.cache/' + cache_key, 'rb') as cached:
                isodist_results = cPickle.load(cached)
        except IOError:
            log_main.info('No pickle cache found. We\'ll write one momentarily.')

    if not isodist_results:
        isodist_results = read_isodist_results('./isodist/', peptide_dict, options.skip_graphs)
        if options.force_cache:
            log_main.debug('Writing isodist result data to pickle cache for future use.')
            with open('.cache/' + cache_key, 'wb') as pickle_cache:
                cPickle.dump(isodist_results, pickle_cache)
    else:
        # We got cached isodist results. Cool!
        log_main.warning('Using cached isodist data: If you\'ve re-run isodist, the cache will be stale!')

    log_main.info("Reading isodist results took: %0.2f secs." % (time.time() - itime))

    # Why not make some matplotlib graphs?
    if options.skip_graphs:
        log_main.warning('Skipping peak graph generation as requested. (Graphs may be stale.)')
    else:
        make_peak_graphs(peptide_dict, isodist_results, results_path)

    # Choose winners: rank predictions and choose the best FRC_NX value
    peptide_predictions, peptide_fail_count, peptide_fail_percent = pick_FRC_NX(peptide_dict, isodist_results)

    del isodist_results  # This is huge.

    # Choose protein-level predictions given the peptides
    protein_predictions, protein_fail_count, protein_fail_percent = pick_protein_enrichment(dta_select_data, peptide_dict,
                                                                                            peptide_predictions)

    # Save output as CSV & HTML.
    generate_output(dta_select_data, peptide_dict,
                    peptide_predictions, peptide_fail_count, peptide_fail_percent,
                    protein_predictions, protein_fail_count, protein_fail_percent,
                    results_path, input_directory)

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
    local_isodist_version = stdout[23:27]  # Snip the isodist version string.

    # We could use pkg_resources.parse_version here, but that isn't always around.
    # Besides, there appears to be only one isodist version for the past ~5 years.
    if local_isodist_version != '2008':
        raise FatalError('isodist is outdated. Please use a version >= 2008')

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
            log_prerun.warning('isodist input file has been modified: %s' % (fname,))

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

    dta_select_csv = csv.reader(open(DTASelect_file, 'rb'), delimiter='\t')

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
                del line[1]  # This is not friendly.

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
    charge_dist = [0] * 6

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

            calc_mh = peptide_data['calc_mh']  # This is the +1 state from the DTASelect file
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
                             if (calc_mz - (MS1_WINDOW / 2.0)) < m < (n15_adjusted_mass + (MS1_WINDOW / 2.0))]

            peptide_dict[peptide_key] = {'sequence': peptide_sequence,
                                         'mz': calc_mz,
                                         'n15mz': n15_adjusted_mass,
                                         'charge': charge,
                                         'peaks': extracted_ms1, }

    extract_peak_logger.info('Summary of charge distribution:')
    for index, item in enumerate(charge_dist):
        extract_peak_logger.info('\t%s: %0.2f%%' % (index, (item / float(sum(charge_dist)) * 100)))

    return peptide_dict


def make_peak_graphs(peptide_dict, isodist_results, results_path):
    """
    Make some graphs of the peaks.

    You could try this with multiprocessing -- I used map_async here for a while, but it seemed to
    unnecessarily complicate things and periodically deadlock. I don't think there's any real
    speedup, given the Python GIL.
    """

    graph_log = logging.getLogger('make_peak_graphs')

    graph_log.debug('Making graph output directories in results/%s/graphs/' % (results_path, ))
    # Let's create output directories for graphs
    for peptide_key in peptide_dict.iterkeys():
        try:
            # Makedirs will make the parent /graphs/ if it doesn't exist.
            os.makedirs('./results/%s/graphs/%s' % (results_path, peptide_key))
        except OSError:
            # Dir might already exist. If it's unwriteable, we'll FATAL it later.
            pass

    graph_log.info('Generating graphs. This will take some time.')
    tasks = [(key, val, isodist_results[key], results_path) for key, val in peptide_dict.iteritems()]
    for key, val in peptide_dict.iteritems():
        graph_log.debug('\tGenerating graph for: %s' % (key,))
        # This shouldn't fail, unless perhaps the output directories aren't writeable?
        _peak_graph_cmd(key, val, isodist_results[key], results_path)

    graph_log.info('Graphs generated successfully.')
    return True


def pick_FRC_NX(peptide_dict, isodist_results):
    """
    Choose the best isodist peptide FRC_NX enrichment guess.
    """

    frc_nx_log = logging.getLogger('pick_FRC_NX')
    frc_nx_log.info('Determining optimal peptide enrichment percentages.')
    peptide_count = len(peptide_dict.keys())
    frc_nx_log.info('Predictions were made for %s peptides.' % (peptide_count,))

    # I'm not sure what the best way to settle on a prediction is. The options include:
    #  - Direct comparsions to raw peaks (ignoring the isodist CHI_SQ scores)
    #  - Averaging/other statistics based on CHI_SQ or other values
    #  - Windowing/mode selection
    #  - K-Means / C-Means / Hierarchical Clustering
    #
    # I've tried a lot of these, and settled on a sorted-windowing approach. This is basically
    # a poor man's divisive hierarchical clustering strategy.
    #
    # The code below does:
    #  1. For the N_PERCENT_RANGE values [0, 10, 20 ...], take each predicted FRC_NX and ...
    #  2. Add it to a list (or heapq if you have a ton of these and want O(lgn) )
    #  3. Loop over the list (pop the heap), looking for values within 1% of each other
    #    - If yes, continue until >= 4 values are within 1% of each other
    #    - If no, move to the next queue entry
    #  4. If >=4 FRC_NX predictions are within 1% of each other, return their mean, otherwise,
    #    we refuse to make an FRC_NX enrichment prediction.
    #
    # WARNING: If you make modifications to N_PERCENT_RANGE, you may wish to tweak this approach.
    peptide_predictions = {}
    fail_count = 0
    for peptide_id in peptide_dict.iterkeys():
        frc_nx_log.debug('Choosing enrichment for %s' % (peptide_id,))
        # We'll fill this dict with each raw prediction plus our
        # computed enrichment value, if we settle on one.
        peptide_predictions[peptide_id] = {}
        isodist_data = isodist_results[peptide_id]

        enrich_list = []
        for percent in N_PERCENT_RANGE:
            # Add raw FRC_NX guesses from isodist
            isodist_guess = isodist_data[percent]['frc_nx']
            peptide_predictions[peptide_id]['guess_' + str(percent)] = isodist_guess
            enrich_list.append(isodist_guess)

        # Print the raw window for debugging.
        frc_nx_log.debug('\tRaw: %s' %
                         ([isodist_data[percent]['frc_nx'] for percent in N_PERCENT_RANGE]))

        predictions_dict = heap_windowing(enrich_list=enrich_list,
                                          margin=0.01,
                                          window_cutoff=4)

        # Did we get any winners? If so, hooray!
        if predictions_dict['golden_window']:
            # Merge in our dict
            for key, value in predictions_dict.iteritems():
                peptide_predictions[peptide_id][key] = value
            frc_nx_log.debug('\tWindow: %s items, %s' % (len(predictions_dict['golden_window']),
                                                         predictions_dict['golden_window']))
            frc_nx_log.debug('\tChoosing: %0.2f%%' % (predictions_dict['guess'] * 100,))
        else:
            frc_nx_log.warning('\tPrediction failed for %s' % (peptide_id,))
            fail_count += 1

    fail_percent = fail_count / peptide_count * 100
    frc_nx_log.info('Prediction failed for %s out of %s peptides (%0.2f%%)' %
                    (fail_count, peptide_count, fail_percent))

    frc_nx_log.info('Peptide enrichment percentages chosen successfully.')
    return (peptide_predictions, fail_count, fail_percent)


def heap_windowing(enrich_list, margin, window_cutoff):
    """
    DESCR GOES HERE
    """
    assert(margin > 0 and margin < 1)
    assert(window_cutoff > 1)

    enrichment_window = []
    golden_window = False
    guess = False

    n = 0
    mean = 0.0
    M2 = 0.0

    for enrich_guess in sorted(enrich_list):
        # This is Welford's algorithm, as implemented by Knuth
        n += 1
        delta = enrich_guess - mean
        mean = mean + (delta / n)
        M2 = M2 + delta * (enrich_guess - mean)

        # Handle the first value
        if len(enrichment_window) == 0:
            enrichment_window.append(enrich_guess)
        # See if we're within MARGIN of the previous value
        elif (enrichment_window[-1] - margin) < enrich_guess < (enrichment_window[-1] + margin):
            enrichment_window.append(enrich_guess)
        # We're not within MARGIN. Clear the window. Let's hope there are more values!
        else:
            enrichment_window = [enrich_guess]

        # Check for (and hang on to) a set of >= window_cutoff enrichment values that are good.
        if len(enrichment_window) >= window_cutoff:
            golden_window = enrichment_window

    # If we got a tiny enrich_list, give up.
    if golden_window and len(enrich_list) > 1:
        guess = sum(golden_window) / float(len(golden_window))

    try:
        variance_n = M2 / n
        variance = M2 / (n - 1)
    except ZeroDivisionError:
        variance = False
        variance_n = False

    return {'golden_window': golden_window,
            'guess': guess,
            'mean': mean,
            'variance': variance,
            'variance_n': variance_n
            }


def pick_protein_enrichment(dta_select_data, peptide_dict, peptide_predictions):
    """
    Pick protein-level enrichments given the peptide enrichment percentages.
    """
    prot_log = logging.getLogger('pick_protein_enrichment')
    prot_log.info('Determining protein-level enrichment percentages.')
    protein_count = len(dta_select_data.keys())
    protein_predictions = {}
    fail_count = 0

    # This is a very similar approach to the above windowed algorithm.
    # Here, we're trying to choose overall protein enrichment given the peptide enrichments.
    for protein_id in dta_select_data.iterkeys():
        enrich_list = []
        # peptide_count = dta_select_data[protein_id]['metadata']['seq_count']
        peptide_count = 0
        for peptide_id in dta_select_data[protein_id]['peptides'].iterkeys():
            peptide_count += 1
            if 'guess' in peptide_predictions[peptide_id]:
                enrich_list.append(peptide_predictions[peptide_id]['guess'])

        # Print the whole list for debugging
        prot_log.debug('\tRaw: %s' % (enrich_list, ))

        predictions_dict = heap_windowing(enrich_list=enrich_list,
                                          margin=0.1,
                                          window_cutoff=2)

        # Make a dict to hold results
        protein_predictions[protein_id] = {}

        # Did we get any winners? If so, hooray!
        if predictions_dict['golden_window']:
            # Merge in our dict
            for key, value in predictions_dict.iteritems():
                protein_predictions[protein_id][key] = value
            protein_predictions[protein_id]['num_samples'] = len(enrich_list)

            prot_log.debug('\tWindow: %s items, %s' % (len(predictions_dict['golden_window']),
                                                       predictions_dict['golden_window']))
            prot_log.debug('\tChoosing: %0.2f%%' % (predictions_dict['guess'] * 100,))
        else:
            prot_log.warning('\tPrediction failed for %s' % (peptide_id,))
            fail_count += 1

    fail_percent = fail_count / protein_count * 100
    prot_log.info('Prediction failed for %s out of %s proteins (%0.2f%%)' %
                  (fail_count, protein_count, fail_percent))

    prot_log.info('Protein enrichment predictions complete.')
    return (protein_predictions, fail_count, fail_percent)


def generate_output(dta_select_data, peptide_dict,
                    peptide_predictions, peptide_fail_count, peptide_fail_percent,
                    protein_predictions, protein_fail_count, protein_fail_percent,
                    results_path, input_directory):
    """
    Save a final output summary of our predictions. This outputs as CSV & HTML files.

    Note: This uses DataTables -- if the number of columns changes, the table may break.
       Be sure to update the HTML headers if you change things.
    """
    output_log = logging.getLogger('generate_output')
    output_log.info('Saving ouput to: results/%s' % (results_path,))

    # Copy over core Bootstrap/Jquery/Datatables files to results dir.
    shutil.copytree('.html/assets/', 'results/%s/assets/' % results_path)

    # Read in our generic HTML and templates.
    with open('.html/header.html') as header_fh, open('.html/footer.html') as footer_fh:
        header_template = string.Template(header_fh.read())
        footer = footer_fh.read()

    #### Peptide raw data output
    output_log.info('Generating peptide output data...')

    # Dictionary keys for our output results.
    peptide_dict_keys = ['sequence', 'charge', 'mz', 'n15mz']
    peptide_predict_keys = ['guess', 'variance', 'variance_n']

   # The "guess_0" -> 100 columns (for a much larger table)
    all_isodist_guesses = ['guess_' + str(percent) for percent in N_PERCENT_RANGE]

    # A function for peptide result printing to reduce repetitive code
    def write_peptide_results(table_head, out_handle, full_details=False):
        # Highlight the top navbar
        if full_details:
            out_handle.write(header_template.safe_substitute(peptide_details='active'))
        else:
            out_handle.write(header_template.safe_substitute(peptide_active='active'))
        # This is a table-specific header, since column counts & titles change.
        out_handle.write(table_head)
        for key in peptide_dict.iterkeys():
            print('<tr><td>' + key + '</td><td>', file=out_handle)
            print('</td><td>'.join([str(peptide_dict[key][x]) for x in peptide_dict_keys]), file=out_handle)
            print('</td><td>', file=out_handle)
            # Make the guesses rounded and prettier
            guess = peptide_predictions[key].get('guess')
            if guess:
                print('</td><td>'.join([str(round(peptide_predictions[key][x] * 100, 2)) for x in peptide_predict_keys]), file=out_handle)
            else:
                print('<i style="color:red">Failed</i></td><td></td><td>', file=out_handle)
            if full_details:
                print('</td><td>', file=out_handle)
                print('</td><td>'.join([str(peptide_predictions[key][x]) for x in all_isodist_guesses]), file=out_handle)
            print('</td></tr>', file=out_handle)
        print('</tbody></table>', file=out_handle)
        out_handle.write(footer)

    # Print the summary table
    with open('.html/peptide_summary_table_head.html') as summary_head:
        peptide_summary_table_head = summary_head.read()
    with open('results/%s/%s' % (results_path, 'peptide_summary.html'), 'wb') as htmlout:
        write_peptide_results(table_head=peptide_summary_table_head,
                              out_handle=htmlout,
                              full_details=False)

    # Print the more detailed table
    with open('.html/peptide_table_head.html') as peptide_head:
        peptide_table_head = peptide_head.read()
    with open('results/%s/%s' % (results_path, 'peptide_details.html'), 'wb') as htmlout:
        write_peptide_results(table_head=peptide_table_head,
                              out_handle=htmlout,
                              full_details=True)

    output_log.debug('Peptide HTML successfully generated.')

    # Add the isodist guesses for our output
    peptide_predict_keys.extend(all_isodist_guesses)

    # Open & write our CSV files
    # Protip: Write to sys.stderr if you're debugging this.
    with open('results/%s/%s' % (results_path, 'all_peptides.csv'), 'wb') as csvout:
        csv_out = csv.DictWriter(csvout, ['pep_key'] + peptide_dict_keys + peptide_predict_keys, extrasaction='ignore')
        csv_out.writeheader()
        for key in peptide_dict.iterkeys():
            csv_out.writerow(dict([('pep_key', key)] + peptide_dict[key].items() + peptide_predictions[key].items()))
    # And TSV, too.
    with open('results/%s/%s' % (results_path, 'all_peptides.tsv'), 'wb') as tsvout:
        tsv_out = csv.DictWriter(tsvout, ['pep_key'] + peptide_dict_keys + peptide_predict_keys, extrasaction='ignore', dialect=csv.excel_tab)
        tsv_out.writeheader()
        for key in peptide_dict.iterkeys():
            tsv_out.writerow(dict([('pep_key', key)] + peptide_dict[key].items() + peptide_predictions[key].items()))
    output_log.debug('Peptide CSV/TSV sucessfully generated.')

    #### End of Peptide data generation

    #### Generate Protein data
    output_log.info('Generating protein output data...')

    prot_prediction_keys = ['guess', 'variance', 'variance_n']
    prot_metadata_keys = ['name', 'validated', 'spect_count', 'molwt', 'length', 'seq_cov', 'pI', 'seq_count']

    # A function for our protein HTML output
    def write_protein_results(table_head, out_handle, full_details=False):
        # Highlight the top navbar
        if full_details:
            out_handle.write(header_template.safe_substitute(protein_details='active'))
        else:
            out_handle.write(header_template.safe_substitute(protein_active='active'))
        # This is a table-specific header, since column counts & titles change.
        out_handle.write(table_head)
        for key in dta_select_data.iterkeys():
            # Without full details, we print protin metadata
            if not full_details:
                print('<tr><td>' + key + '</td><td>', file=out_handle)  # Protein Key
                print('</td><td>'.join([str(dta_select_data[key]['metadata'][x]) for x in prot_metadata_keys]), file=out_handle)
                print('</td><td>', file=out_handle)
                # Beautify the output
                guess = protein_predictions[key].get('guess')
                if guess:
                    print('</td><td>'.join([str(round(protein_predictions[key][x] * 100, 2)) for x in prot_prediction_keys]), file=out_handle)
                else:
                    print('<i style="color:red">Failed</i></td><td></td><td>', file=out_handle)
            else:
                # Print all individual peptides.
                print('</td></tr>', file=out_handle)
                peptides_from_key = dta_select_data[key]['peptides']
                for peptide_id in peptides_from_key.iterkeys():
                    print('<tr><td>' + key + '</td><td>' + peptides_from_key[peptide_id]['sequence'] + '</td><td>', file=out_handle)
                    peptide_guess = peptide_predictions[peptide_id].get('guess')
                    if peptide_guess:
                        print('</td><td>'.join([str(round(peptide_predictions[peptide_id][x] * 100, 2)) for x in prot_prediction_keys]), file=out_handle)
                    else:
                        print('<i style="color:red">Failed</i></td><td></td><td>', file=out_handle)
            print('</td></tr>', file=out_handle)
        print('</tbody></table>', file=out_handle)
        out_handle.write(footer)

    # Write the protein-level summary page
    with open('.html/protein_summary_table_head.html') as protein_head:
        protein_summary_table_head = protein_head.read()
    with open('results/%s/%s' % (results_path, 'protein_summary.html'), 'w') as by_protein:
        write_protein_results(table_head=protein_summary_table_head,
                              out_handle=by_protein,
                              full_details=False)

    # And the more detailed page
    with open('.html/protein_table_head.html') as protein_head:
        protein_table_head = protein_head.read()
    with open('results/%s/%s' % (results_path, 'protein_details.html'), 'w') as by_protein:
        write_protein_results(table_head=protein_table_head,
                              out_handle=by_protein,
                              full_details=True)

    output_log.debug('Protein HTML successfully generated.')

    # Generate CSV/TSV output
    prot_output_keys = ['sequence', 'charge', 'mz', 'n15mz', 'guess', 'guess']

    with open('results/%s/%s' % (results_path, 'all_proteins.csv'), 'wb') as csvout:
        csv_out = csv.DictWriter(csvout, ['prot_key'] + prot_metadata_keys + prot_prediction_keys, extrasaction='ignore')
        csv_out.writeheader()
        for key in dta_select_data.iterkeys():
            csv_out.writerow(dict([('prot_key', key)] + dta_select_data[key]['metadata'].items() + protein_predictions[key].items()))
    # And TSV, too.
    with open('results/%s/%s' % (results_path, 'all_proteins.tsv'), 'wb') as tsvout:
        tsv_out = csv.DictWriter(tsvout, ['prot_key'] + prot_metadata_keys + prot_prediction_keys, extrasaction='ignore', dialect=csv.excel_tab)
        tsv_out.writeheader()
        for key in dta_select_data.iterkeys():
            tsv_out.writerow(dict([('prot_key', key)] + dta_select_data[key]['metadata'].items() + protein_predictions[key].items()))

    output_log.debug('Protein CSV/TSV sucessfully generated.')

    output_log.info('Finalizing HTML output statistics...')
    # And the index template:
    with open('.html/index_templ.html') as index_fh:
        index_template = string.Template(index_fh.read())

    # Metadata for the index page.
    index_template_keys = {'run_id': results_path,
                           'run_date': datetime.datetime.now().strftime("%m/%d/%Y"),
                           'input_directory': input_directory,
                           'protein_count': len(dta_select_data.keys()),
                           'protein_success': len(dta_select_data.keys()) - protein_fail_count,
                           'protein_percent': 100 - protein_fail_percent,
                           'peptide_count': len(peptide_dict.keys()),
                           'peptide_success': len(peptide_dict.keys()) - peptide_fail_count,
                           'peptide_percent': 100 - peptide_fail_percent,
                           }
    with open('results/%s/%s' % (results_path, 'index.html'), 'w') as index:
        index.write(header_template.safe_substitute(index_active='active'))
        index.write(index_template.substitute(index_template_keys))
        index.write(footer)

    output_log.info('Results successfully written.')
    return True


def _peak_graph_cmd(peptide_key, peptide_value, isodist_results, results_path):
    """
    Generate peak graphs from make_peak_graphs.
    """
    # Not that this app is secure, but just so we don't break things, sanitize some outputs.
    valid_chars = "-_.()%s%s" % (string.ascii_letters, string.digits)

    # Semi-sanitized. Not secure, but we aren't likely to accidentally break things if peptide_key is insane.
    safer_key = ''.join(c for c in peptide_key if c in valid_chars)

    # Our m/z peaks don't change through the loop. Just the predictions.
    m, z = zip(*peptide_value['peaks'])

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
        ax.plot(m, z, 'r-', linewidth=0.5)

        # Get the m/z data.
        # Note: These are pre-filtered to exclude all pairs with intensity <=0, as they mess w/ the graph.
        isodist_m = [m / float(peptide_value['charge']) for m in isodist_data['peak_fit'][0]]
        isodist_z = isodist_data['peak_fit'][1]

        ax.plot(isodist_m, isodist_z, 'b-.', linewidth=1.2)
        #ax.autoscale_view() # I really don't know what this does.
        ax.grid(True)

        matplotlib.pyplot.savefig("results/%s/graphs/%s/%s.png" % (results_path, safer_key, n_percent))

        # I'm not sure if this matters for memory usage.
        # Matplotlib uses a *lot* of RAM and nothing seems to reduce it.
        fig.clf()
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
                        'n_percent': 0000, }  # Note: n_percent is adjusted on-the-fly below.

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
        isodist_log.info('Running isodist jobs locally.')

        pool = multiprocessing.Pool(num_threads)
        #tasks = peptide_dict.keys()
        tasks = [x + "/" + str(y) for x in peptide_dict.keys() for y in N_PERCENT_RANGE]
        results = []

        # TODO: Use error_callback?)
        r = pool.map_async(_isodist_cmd, tasks, num_threads, callback=results.append)
        # We use .get() instead of .wait() due to some interrupt issues (Fixed in Py3.3?)
        # This doesn't totally fix the proble, but it makes things a little better.
        # See:
        #  http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
        #  http://bugs.python.org/issue8296
        r.get(999999)  # Block until our pool returns
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
    read_logger.info('Reading isodist results. This will be slow if you are generating graphs.')

    # Save our results in a dict
    isodist_results = {}

    # We'll save some of the .batch.csv results
    keys = ['mw', 'z', 'tim', 'chisq', 'symb', 'mz', 'b', 'off', 'gw', 'amp_u', 'amp_l', 'frc_nx']
    casts = [float, int, str, float, int, float, float, float, float, float, float, float]
    # We ran isodist on each peptide

    read_logger.debug('Reading %s percent values over %s keys.' % (len(N_PERCENT_RANGE), len(peptide_dict)))
    progress = 0
    total = len(N_PERCENT_RANGE) * len(peptide_dict)

    for peptide_key, n_percent in [(x, y) for x in peptide_dict.iterkeys() for y in N_PERCENT_RANGE]:
        progress += 1
        read_logger.debug('%0.2f%% complete. On key: %s, %s%%' % ((progress / total) * 100, peptide_key, n_percent))

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
                peak_fit = csv.reader(fit_file, quoting=csv.QUOTE_NONNUMERIC)  # Pre-cast to float
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
        stdoutdata, stderrdata = p.communicate()  # This will block until isodist finishes.
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
    parsed_xml = cElementTree.parse(mzXML_file)
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

            decoded_b64 = b64decode(scan[0].text)  # Decode the packed b64 raw peaks
            # These are packed as big-endian IEEE 754 binary32
            floating_tuple = unpack('>' + str(len(decoded_b64) // 4) + 'f', decoded_b64)
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
