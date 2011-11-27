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
from ruffus import *
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

# These paths point to the Bowtie indices for these reference genomes.
# They are used by filter_human_contaminant and filter_mouse_contaminant
HG19_BOWTIE2_INDEX = "/home/comp/jglab/semenko/databases/bowtie-indexes/hg19/hg19"
MM9_BOWTIE2_INDEX = "/home/comp/jglab/semenko/databases/bowtie-indexes/mm9/mm9"

# Path to the "BigBelly" database of ~all gut sequences.
# This contains sequences merged from NCBI-gut NCBI-gut-DRAFT and HMRGD
BIG_BELLY_BOWTIE2_INDEX = ""
# Path to directory of annotations (KEGG/COG, etc.) of BigBelly DB.
BIG_BELLY_ANNOTATIONS_DIR = ""

## *****
## Config Options
## *****

# Should we submit jobs to DRMAA? If false, we just run locally. This depends on DRMAA.
USE_DRMAA = True

# If a program has thread options (e.g. bowtie2 with --threads <int>), how many
# threads should we spawn?
# Of note -- Older versions of SGE may have issues with allocating >1 core/job.
MAX_THREADS = 1

# When should we warn a user about a possibly "contaminated" sample (with HG19/MM9)?
# These default settings are arbitrary. A contaminant flag will only warn the user -- the pipeline will still run.
CONTAMINANT_PERCENT_TRIGGER = 5.0
CONTAMINANT_READ_COUNT_TRIGGER = 15000

# These set paths & minimum version requirements. Checked in TODO
# Keeping a tool somewhere else? Try: 'bowtie2': ('/my/path/to/bowtie2', '2.0.0-beta3')
TOOLS_AND_VERSIONS = {
    'bowtie2': ('bowtie2', '2.0.0-beta3'),
    'NCBI BLAST+': ('blastn', '2.2.25'),
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
    """The main part of our metagenomic pipeline."""
    global USE_DRMAA, DRMAA_SESSION, CMD_LINE_INPUT_FILES

    # Begin timing execution

    starttime = time.time()

    parser = OptionParser(usage="usage: %prog [options] input_file(s)\n\nInput files may be fasta, fastq, or SFF.",
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

    if not len(args):
        parser.error("You must specify one or more input files (.sff, .fa, .fq, ...).\n\nTry --help for help.")

    CMD_LINE_INPUT_FILES = args

    for input_file in CMD_LINE_INPUT_FILES:
        if not os.access(input_file, os.R_OK):
            parser.error("Cannot read input file.")
        # TODO: Remove this, or auto-detect filetype (using mimetypes)?
        if not input_file.endswith(('.fa', '.fna', '.ffn', '.fasta', '.fastq', '.fq', '.sff')):
            parser.error("Input filetype not recognized (make sure your files end in something sensible, like .fastq).")

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

    # Let's cook this turkey!
    pipeline_run([generate_run_summary], verbose = 5)

    # Cleanup.
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
    # Make sure bowtie2 exists
    if not len(which(TOOLS_AND_VERSIONS['bowtie2'][0])):
        raise FatalError('bowtie2 not found or not executable.')

    # Check bowtie2 version
    bowtie_version_process = subprocess.Popen([which(TOOLS_AND_VERSIONS['bowtie2'][0]), '--version'],
                                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bowtie_version_process.wait()
    _, stderr = bowtie_version_process.communicate()
    local_bowtie_version = stderr[17:28]

    if not cmp(parse_version(local_bowtie_version), parse_version(TOOLS_AND_VERSIONS['bowtie2'][1])) >= 0:
        raise FatalError('bowtie2 is outdated. Please use a version >= %s' % TOOLS_AND_VERSIONS['bowtie2'][1])

    # bowtie2 database checks
    if not os.access(HG19_BOWTIE2_INDEX + ".1.bt2", os.R_OK):
        raise FatalError('Unable to read HG19 bowtie2 index. Have you edited the config options in this script?')
    if not os.access(MM9_BOWTIE2_INDEX + ".1.bt2", os.R_OK):
        raise FatalError('Unable to read MM9 bowtie2 index. Have you edited the config options in this script?')

    return True


def _input_file_decorator_helper():
    """
    Because decorators are interpreted early-on, we use this function to pass the input
    to Ruffus' @files decorator. See: http://www.ruffus.org.uk/decorators/files_ex.html
    """
    global CMD_LINE_INPUT_FILES

    # TODO: Clean this up.
    # TODO: Make sure we check for collisions early (e.g. abc.fa and abc.fasta both passed.)

    # Make a directory for the basename of each input file
    # e.g. "my.sequence.123.fastq" -> "my.sequence.123/"
    for element in CMD_LINE_INPUT_FILES:
        try:
            os.mkdir(os.path.splitext(element)[0])
        except OSError:
            # Directory already exists.
            pass

    return [[element, os.path.splitext(element)[0] + '/' + element + '.qc-filtered.gz'] for element in CMD_LINE_INPUT_FILES]

###########################
# Actual pipeline code is here.
# DO NOT MODIFY unless you understand Python decorators.
#
# This code depends on Ruffus (MIT License), see http://code.google.com/p/ruffus/
##########################

@files(_input_file_decorator_helper)
def quality_control_filter(input_file, output_file):
    """
    QC filtering on our input files.
    """
    qc_log = logging.getLogger('quality_control_filter')
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


# Parallel contaminant filtering [Human, Mouse, & pathogen detection]
@follows(quality_control_filter)
@transform(quality_control_filter, suffix('.qc-filtered.gz'), '.human-filtered-only.gz')
def filter_human_contaminant(input_file, output_file):
    pass

@follows(quality_control_filter)
def filter_mouse_contaminant(input_file, output_file):
    pass

@follows(quality_control_filter)
def stringent_pathogen_detection(input_file, output_file):
    pass



# Linear block -- this short segment is the rate limiting step
@follows(filter_human_contaminant, filter_mouse_contaminant, stringent_pathogen_detection)
def duplicate_detection(input_file, output_file):
    pass

@follows(duplicate_detection)
def align_to_gut_database(input_file, output_file):
    pass



# Ok, back to parallel tasks!
# Let's get the interesting database reads from the aligned sequences (KEGG, COG, etc.)
@follows(duplicate_detection)
def generate_taxonomy(input_file, output_file):
    pass


# Let's also take the un-aligned reads and see what we can do with them.
@follows(align_to_gut_database)
def align_leftovers_to_viral_db():
    pass

@follows(align_to_gut_database)
def align_leftovers_to_kegg():
    pass

@follows(align_to_gut_database)
def align_leftovers_to_cog():
    pass




# Create a final run summary and declare victory!
@follows(align_leftovers_to_viral_db, align_leftovers_to_kegg, align_leftovers_to_cog)
def generate_run_summary():
    pass


# Generate a dependency graph, if asked
def make_dependency_graph():
    pipeline_printout_graph("flowchart.svg", "svg", [generate_run_summary], no_key_legend = True)

#make_dependency_graph()


### ---------------------------------------------
### Helper Functions
### ---------------------------------------------

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

    # Per Ruffus FAQ, we have a random delay so SGE has a chance get its emotions in check.
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
