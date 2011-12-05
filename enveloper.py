#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright (c) 2011 Nick Semenkovich <semenko@alum.mit.edu> / WUSTL
# 
# Developed for the Gordon Lab, Washington University in St. Louis (WUSTL)
# http://gordonlab.wustl.edu/
#
# This software is released under the MIT License <http://www.opensource.org/licenses/mit-license.php>
#
"""
Some envelope stuff with mass spec stuff.
"""

__author__ = 'Nick Semenkovich <semenko@alum.mit.edu>'
__copyright__ = 'Gordon Lab at Washington University in St. Louis / gordonlab.wustl.edu'
__license__ = 'MIT'
__version__ = '1.0'

from optparse import OptionParser, OptionGroup
from pkg_resources import parse_version, get_distribution
import csv
import datetime
import socket
import logging
import os
import subprocess
import random
import time
import signal
import sys

### ---------------------------------------------
### Configuration Variables.
### You may need to edit this block for your setup.
### ---------------------------------------------


## *****
## Data Paths
## *****



## *****
## Config Options
## *****

# Should we submit jobs to DRMAA? If false, we just run locally. This depends on DRMAA.
USE_DRMAA = True

# If a program has thread options (e.g. bowtie2 with --threads <int>), how many
# threads should we spawn?
# Of note -- Older versions of SGE may have issues with allocating >1 core/job.
MAX_THREADS = 1

# These set paths & minimum version requirements. Checked in TODO
# Keeping a tool somewhere else? Try: 'bowtie2': ('/my/path/to/bowtie2', '2.0.0-beta3')
TOOLS_AND_VERSIONS = {
    'isodist': ('isodist', '2008'),
    'other': ('tool', '1.0'),
}

### ---------------------------------------------
### Exceptions & Handlers
### ---------------------------------------------

def handle_SIGINT(signal, frame):
    raise FatalError('User abort: Caught SIGINT')

class FatalError(Exception):
    """ Thrown when we should die. """
    global USE_DRMAA, DRMAA_RUNNING_JOBS

    def __init__(self, message):
        Exception.__init__(self)

        fatal_logger = logging.getLogger('FatalError')

        if USE_DRMAA:
            # Terminate any running jobs.
            for job_id in DRMAA_RUNNING_JOBS:
                DRMAA_SESSION.control(job_id, drmaa.JobControlAction.TERMINATE)
                fatal_logger.critical('Killed SGE/DRMAA job %s' % job_id)

        fatal_logger.error(message)
        print "Terminating ..."
        sys.exit(1)

### ---------------------------------------------
### Colorize the logger class
### Via StackOverflow: http://stackoverflow.com/questions/384076/how-can-i-make-the-python-logging-output-to-be-colored
### ---------------------------------------------
class ColorFormatter(logging.Formatter):
  FORMAT = ("$BOLD%(name)-25s$RESET: %(levelname)-20s: "
            "%(message)-100s "
            "($BOLD%(filename)s$RESET:%(lineno)d)")

  BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)

  RESET_SEQ = "\033[0m"
  COLOR_SEQ = "\033[1;%dm"
  BOLD_SEQ = "\033[1m"

  COLORS = {
    'WARNING': YELLOW,
    'INFO': WHITE,
    'DEBUG': BLUE,
    'CRITICAL': YELLOW,
    'ERROR': RED
  }

  def formatter_msg(self, msg, use_color = True):
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

    # OPT:
    # file hint (sff/fa/fq)
    # split on MIDs.

    # This ignores hyperthreading pseudo-cores, which is fine since we presumably hose the ALU.
    # TODO: Change default when DRMAA is enabled.
    parser.add_option("--threads", help="Threads. We parallelize the invocations of RNAhybrid. [Default: # of CPU cores]",
                      default=os.sysconf('SC_NPROCESSORS_ONLN'), action="store", type="int", dest="threads")
    parser.add_option("--profile", help="Profile. Invokes the yappi python profiling engine. Will slow down execution.",
                      default=False, action="store_true", dest="profileMe")

    group = OptionGroup(parser, "Range Settings (optional)")
    group.add_option("--start-num", help="What number miRNA ortholog group to start scanning from (inclusive).",
                      default=-1, action="store", type="int", dest="startNum")
    group.add_option("--stop-num", help="What number miRNA ortholog group to STOP scanning at (exclusive).",
                      default=-1, action="store", type="int", dest="stopNum")
    parser.add_option_group(group)

    
    (options, args) = parser.parse_args()

    # Sanity check range inputs
    range_given = False
    if options.startNum >= 0 or options.stopNum >= 0:
        if not (options.startNum >= 0 <= options.stopNum):
            parser.error("If you specifiy a start/stop, you must specify both ends of the range!")
        if options.startNum >= options.stopNum:
            parser.error("Invalid scan range.")
        if options.startNum == 1:
            print "WARNING: This is a Python range, where lists are zero indexed! Are you sure you mean '1'?"
        range_given = True

    if not len(args) == 1:
        parser.error("You must specify exactly one input directory.\n\nTry --help for help.")

    input_directory = args[0]

    if not os.access(input_directory, os.R_OK):
        parser.error("Cannot read input directory.")

    # Make sure we have some appropriately named input files, so we don't die later.
    directory_list = os.listdir(input_directory)
    if 'DTASelect-filter.txt' not in directory_list:
        parser.error("DTASelect-filter.txt not found in input directory.")
    if not len([file for file in directory_list if file.endswith('.mzXML')]) == 1:
        parser.error("Exactly one .mzXML file must be present in the input directory.")
    if not len([file for file in directory_list if file.endswith('.sqt')]) >= 1:
        parser.error(".sqt file(s) not found in input directory.")


    # Let's set up a logging system
    # We log DEBUG and higher to log file, and write INFO and higher to console.
    datestamp = datetime.datetime.now().strftime("%m%d-%H%M")
    logfile = datestamp + '.' + socket.gethostname().split('.')[0] + '.log'
    logging.basicConfig(filename = logfile, filemode = 'w',
                        format = '%(asctime)s: %(name)-25s: %(levelname)-8s: %(message)s',
                        level = logging.DEBUG)

    # Define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # Point to our colorized formatter.
    console.setFormatter(ColorFormatter())
    logging.getLogger('').addHandler(console)

    log_main = logging.getLogger('main')
    
    log_main.info('Logging to %s' % logfile)

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

    # Take the unique DTA Select peptide fragments and run isodist
    #run_isodist(input_directory.rstrip('/'))

    # Now that we have the isodist predicted spectra, parse the mzXML
    #mzXML_data = parse_mzXML()

    # Let's cook this turkey!
    # (actually do comparisons from isodist <-> mzXML spectra)
    #turkey_cooker(mzXML_data, input_directory.rstrip('/'))

    #pipeline_run([generate_run_summary], verbose = 5)

    # Cleanup.
    if USE_DRMAA:
        DRMAA_SESSION.exit()


def pre_run_version_checks():
    """
    Check that our dependencies are available and meet minimum version requirements.
    """
    log_prerun = logging.getLogger('pre_run_version_checks')

    if sys.hexversion < 0x02060500:
        raise FatalError('Outdated Python version. Please use >=2.6.5')

    if not cmp(parse_version(get_distribution("ruffus").version), parse_version('2.2')) >= 0:
        raise FatalError('Outdated Ruffus version. Please use >= 2.2')

    ### Tool Checks
    # Make sure isodist exists
    if not len(which(TOOLS_AND_VERSIONS['isodist'][0])):
        raise FatalError('isodist not found or not executable.')

    # Check isodist version (only one 2008 version at time of writing)
    isodist_version_process = subprocess.Popen([which(TOOLS_AND_VERSIONS['isodist'][0])],
                                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    isodist_version_process.wait()
    stdout, _ = isodist_version_process.communicate()
    local_isodist_version = stdout[23:27] # Snip the isodist version string.

    if not cmp(parse_version(local_isodist_version), parse_version(TOOLS_AND_VERSIONS['isodist'][1])) >= 0:
        raise FatalError('isodist is outdated. Please use a version >= %s' % TOOLS_AND_VERSIONS['isodist'][1])

#    if not os.access(MM9_BOWTIE2_INDEX + ".1.bt2", os.R_OK):
#        raise FatalError('Unable to read MM9 bowtie2 index. Have you edited the config options in this script?')

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
    peptide_list = []

    # These files aren't the easiest to parse. We pick lines based on the number of TSV elements.
    for line in dta_select_csv:
        print line
        # Have we gone past the header in the file yet?
        if past_header:
            if len(line) == 9:
                # Length 9 lines are protein groups (they delimit peptide sections)
                if added_peptides == False:
                    # This section has multiple protein headers, or we've just started to parse the file.
                    current_keys.append(line[0])
                else:
                    # We must've just entered a new protein section. Reset our "current_keys" list.
                    added_peptides = False
                    current_keys = [line[0]]
                    peptide_list = []

                if line[0] in dta_dict:
                    # I don't think this is ever possible. But let's be paranoid.
                    raise FatalError('Duplicate protein key in DTA Select file!')

                # Make a dict in our dict!
                dta_dict[line[0]] = {'metadata': {}, 'peptides': []}

                # Set our metadata. In the file, the format is:
                # 'Sequence Count', 'Spectrum Count', 'Sequence Coverage', 'Length', 'MolWt', 'pI', 'Validation Status', 'Descriptive Name'
                metadata_dict = dta_dict[line[0]]['metadata']
                keys = ['seq_count', 'spect_count', 'seq_cov', 'length', 'molwt', 'pI', 'validated', 'name']
                type = [int, int, str, int, int, float, str, str]
                for key, value, cast in zip(keys, line[1:], type):
                    metadata_dict[key] = cast(value)
 

            elif len(line) == 13:
                # Length 13 lines are peptide data, which we add to a protein entry.
                
                # Mark that we've started adding data. Otherwise we'd add it to the wrong key!
                added_peptides = True

                # In the file, the format is:
                # ['Unique', 'FileName', 'XCorr', 'DeltCN', 'Conf%', 'M+H+', 'CalcM+H+', 'TotalIntensity', 'SpR', 'Prob Score', 'IonProportion', 'Redundancy', 'Sequence']
                # TODO: Make this suck less.
                for key in current_keys:
                    dta_dict[key]['peptides'].append(line)

            elif len(line) == 4:
                # We're at the end of the file. Victory is ours!
                break
            else:
                raise FatalError('Odd structure in DTA Select file!')
                
        # We aren't entirely though the header yet.
        else:
            # The header is 13 elements and starts with "Unique"
            if len(line) == 13 and line[0] == "Unique":
                past_header = True
                parse_dta_log.info('Found data stream in DTA Select file.')


    parse_dta_log.info('Finished parsing.')

    print dta_dict

    return dta_dict


def run_isodist(output_path, isotope_distributions):
    """
    Run isodist.
    """
    run_isodist_log = logging.getLogger('run_isodist')
    run_isodist_log.info('Running isodist.')
    # Should we spawn jobs via DRMAA?
    if USE_DRMAA == True:
        pass
    else:
        pass


def parse_mzXML(mzXML_file):
    """
    Open and parse an mzXML file. TODO: Object orient this or something.
    """
    parse_mzXML_log = logging.getLogger('parse_mzXML')
    # TODO: Run FastQC?
    # fasta_formatter -i 50k.fna -w 0 | fastx_trimmer -m 60 -t 60
    # We accept fasta/fastq/sff, so it'd better be one of those formats.
    if input_file.endswith(('.fa', '.fna', '.ffn', '.fasta')):
        # FASTA input
        qc_log.info('QC filtering FASTA file %s (Output: %s)' % (input_file, output_file))
        return_code = subprocess.call("fasta_formatter -i 50k.fna -w 0 | fastx_trimmer -m 60 -t 60", shell=True)
        if return_code != 0:
            raise FatalError('QC Filtering failed.')
        
    elif input_file.endswith(('.fastq', '.fq')):
        # FastQ
        qc_log.info('QC filtering FASTQ file %s (Output: %s)' % (input_file, output_file))
    elif input_file.endswith('.sff'):
        # SFF
        qc_log.info('QC filtering SFF file %s (Output: %s)' % (input_file, output_file))
    else:
        # We should never get here, since main() checks input extensions.
        raise FatalError('Unknown input file extension.')


# Very annoying that Python doesn't have a `which` equivalent. I avoid calling sys(which), since,
# well, you /could/ run this on Windows. :/
# Adapted from: http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
def which(program):
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def deploy_drmaa_job(job_command, job_parameters):
    global DRMAA_RUNNING_JOBS

    drmaa_logger = logging.getLogger('deploy_drmaa_job')

    assert(os.access(job_command, ox.X_OK))
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

    # Cleanup.
    # TODO: Use runBulkJobs when parallelizing?
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
