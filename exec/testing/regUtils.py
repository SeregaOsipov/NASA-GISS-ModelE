# This module contains utilities to help setup the regression tests
import string
import ConfigParser
import os
import re
import sys
import errno
import datetime
import shutil
import subprocess
import logging
import regTest

logger = logging.getLogger('utils')

#-------------------------------------------------------------------------------
# Read a config file, return config object
def readConfig(cfgfile):
    config = ConfigParser.ConfigParser()
    config.read(cfgfile)

    logger.info('Read configuration file %s',cfgfile)
    return config

#-------------------------------------------------------------------------------
# Create/return a list of model run configurations specified in config file
def getModelConfigurations(config):
    sections = config.sections()
    modelConfig = {}

# Retrieve all model run configurations. For convenience divide the sections
# in the configuration file into two types: CONFIG and others. The former 
# have CONFIG in their names. Thus, if a section name does NOT have CONFIG 
# in its name then it is a model configuration.
    for sect in sections:
        match = not re.search("CONFIG",sect)
        # get all model config sections from config file
        if (match):
            modelConfig[sect] = ConfigSectionMap(config, sect)

# Store model configurations in a list and let each item in the list
# have access to the user-defined options
    runList = []
    for name,options in modelConfig.items():
        # Each item is a regression test (regTest) instance
        runList.append(regTest.regTest(name))

    for d in runList:
        d.setOpts(modelConfig)
        
    return runList

#-------------------------------------------------------------------------------
# Return a dict (i.e. a key,value pair) from each section in config
def ConfigSectionMap(config, section):
    adict = {}
    options = config.options(section)

    for option in options:
        try:
            adict[option] = config.get(section, option)
            if adict[option] == -1:
                logger.info('skip option %s',option)
        except:
            logger.error("exception on %s!" % option)
            adict[option] = None

    logger.debug('Read section %s',section)

    return adict

#-------------------------------------------------------------------------------
# Get a list of compilers used
def getCompilers(config):
   compconfig = ConfigSectionMap(config, 'COMPCONFIG')
   complist = compconfig['compilers'].split(",")
   compilers = []
   for c in complist:
       compilers.append(c.strip())
   return compilers

#-------------------------------------------------------------------------------
# Return a command that creates a clone of the reference clone
def gitCloneCommand(config, expname, compiler, cmode):
    userconfig = ConfigSectionMap(config, 'USERCONFIG')
    branch    = userconfig['repobranch']
    if not branch:
        branch = 'detached'
    scratch   = userconfig['scratchdir'] + '/scratch/' + branch
    reference = scratch + '/' + branch
    clone     = scratch + '/' + compiler + '/' + expname + '.' + cmode

    if not os.path.isdir(clone):
        if branch == 'detached':
            s = string.Template('git clone $r $t > /dev/null 2>&1')
            return s.substitute(r=reference, t=clone)
        else:
            s = string.Template('git clone -b $b $r $t > /dev/null 2>&1')
            return s.substitute(b=branch, r=reference, t=clone)
    else:
        logger.debug('Git clone %s exists', clone)
        return clone

#-------------------------------------------------------------------------------
# Workaround for "mkdir -p" command
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

#-------------------------------------------------------------------------------
# "Safe" way to clean the contents of a directory
def cleanDir(adir):
    logger.info('Cleaning up scratch space...')
    if(adir == '/' or adir == "\\"):
        logger.error('Cannot clean %s',adir)
        return
    else:
        for file_object in os.listdir(adir):
            logger.debug('Will clean up %s',adir)
            file_object_path = os.path.join(adir, file_object)
            if os.path.isfile(file_object_path):
                os.unlink(file_object_path)
            else:
                shutil.rmtree(file_object_path)


#-------------------------------------------------------------------------------
#  Test if an executable program exists in the path - like unix's which
def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
   
#-------------------------------------------------------------------------------
# Create directory with timestamp
def mkdirTimeSTamp(list, filename):
    mydir = os.path.join(os.getcwd(), datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S'))
    try:
        os.makedirs(mydir)
    except OSError, e:
        if e.errno != 17:
            raise # This was not a "directory exist" error..
    with open(os.path.join(mydir, filename), 'w') as d:
        d.writelines(list)

#-------------------------------------------------------------------------------
def copyFile(src, dest):
    try:
        shutil.copy(src, dest)
        return 0
    # eg. src and dest are the same file
    except shutil.Error as e:
        print('Error: %s' % e)
    # eg. source or destination doesn't exist
    except IOError as e:
        print('Error: %s' % e.strerror)
    # eg. Operation not permitted
    except OSError as e:
        print('Error: %s' % e.strerror)
    return 1

#-------------------------------------------------------------------------------
def checkpointName(name, mode, endTime, npes):
#  Return a checkpoint file name with various identifiers
    if mode == 'serial':
        return name  + '.' + endTime
    else:
        return name + '.' + endTime + '.np=' + str(npes)

#-------------------------------------------------------------------------------
def loggerSetup(cfg):
# Logger setup       
    logging.basicConfig(
        filename = str(cfg) + '.LOG',
        format = "%(levelname) -10s %(module)s:%(lineno)s %(funcName)s %(message)s",
        level = logging.DEBUG,
        filemode = 'w'
    )
    stdoutLog = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(name)s : %(message)s')
    stdoutLog.setFormatter(formatter)
    if os.environ.has_key('DEBUG'):
        stdoutLog.setLevel(logging.DEBUG)
    else:
        stdoutLog.setLevel(logging.INFO)
    logger = logging.getLogger()
    logger.addHandler(stdoutLog)

#-------------------------------------------------------------------------------
from ConfigParser import SafeConfigParser
def showConfig(cfgfile):

    parser = SafeConfigParser()
    parser.read(cfgfile)

    for section_name in parser.sections():
        print 'Section:', section_name
        print '  Options:', parser.options(section_name)
        for name, value in parser.items(section_name):
            print '  %s = %s' % (name, value)
        print

#-------------------------------------------------------------------------------
def cleanScratch(config):
    logger.info('Clean up testing environment')
    userconfig = ConfigSectionMap(config, 'USERCONFIG')
    branch    = userconfig['repobranch']
    if not branch:
        branch = 'detached'
    resultsDir = userconfig['scratchdir'] + '/results/' + branch
    scratchDir = userconfig['scratchdir'] + '/scratch/' + userconfig['repobranch']

    if not os.path.exists(resultsDir):
        mkdir_p(resultsDir)    
        mkdir_p(scratchDir)
    else:
        cleanDir(scratchDir)
        cleanDir(resultsDir)

#-------------------------------------------------------------------------------
def header(cfg):
    logger.info('USER CONFIGURATION:')
    userconfig  = ConfigSectionMap(cfg, 'USERCONFIG')
    mailto     = userconfig['mailto']
    branch     = userconfig['repobranch']
    if not branch:
        branch = 'detached'
    buildtype  = userconfig['buildtype']
    sections = cfg.sections()

    print('-'*80)
    print('Testing repository: ' + userconfig['repository'])
    print('Testing _'+branch+'_ branch using _'+buildtype+'_ build type')
    ## To print w/o a CR  append a comma after the last argument to print. 
    ##print 'Rundecks: ',
    ##for section_name in sections:
    ##    if 'CONFIG' not in section_name:
    ##        print section_name,
    ##print
    print('Output in: ' +  userconfig['scratchdir'])
    if mailto:
        print('Results will be mailed to: ' + mailto)
    print('-'*80)

#-------------------------------------------------------------------------------
def writeDiff(results, fileH):
    fileH.write('%20s' % (results[0]))
    fileH.write('%10s' % (results[1]))
    fileH.write('%8s'  % (results[2]))
    fileH.write('%4s'  % '    ')
    for s in results[3:]:
        fileH.write('{: ^5}'.format(s))
        fileH.write('%3s'  % '   ')
    fileH.write('\n')

      

