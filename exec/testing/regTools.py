# This module contains tools to help setup the modelE regression tests
import string
import ConfigParser
import re
import os
import sys
import shutil
import subprocess as sp
import glob
import logging
import time
import regUtils as util

logger = logging.getLogger('tools')

#-------------------------------------------------------------------------------
# Setup modelE testing environment:
def setupEnv(config, compconfig):
    logger.info('Setup testing environment')
    userconfig = util.ConfigSectionMap(config, 'USERCONFIG')
    branch =  userconfig['repobranch']
    if not branch:
        branch = 'detached'
    resultsDir = userconfig['scratchdir'] + '/results/' + branch
    scratchDir = userconfig['scratchdir'] + '/scratch/' + branch
    makesystem =  userconfig['makesystem']

    # Make sure - if specified - that work space is clean
    if userconfig['cleanscratch'] == 'yes':
        util.cleanScratch(config)

        if makesystem == "makeOld":
           setupModelEenv(config, compconfig)
        gitCloneRepository(config)


#-------------------------------------------------------------------------------
# Clone the model from the user-specified git repository
def gitCloneRepository(config):
    logger.info('Clone user-specified git repository')
    userconfig = util.ConfigSectionMap(config, 'USERCONFIG')
    scratch = userconfig['scratchdir']
    repo = userconfig['repository']
    branch =  userconfig['repobranch']

    if not branch:
        branch = 'detached'
        clone = scratch + '/scratch/' + branch + '/' + branch
        cmd = (['git', 'clone', repo, clone])
    else:
        clone = scratch + '/scratch/' + branch + '/' + branch
        cmd = (['git', 'clone', '-b', branch, repo, clone])

    resultsDir = userconfig['scratchdir'] + '/results/' + branch + '/'

    cwd = os.getcwd()
    logger.debug('Cloning %s into %s', repo, clone)

    proc = sp.Popen(cmd)
    proc.wait()
	
    os.chdir(clone)
    cmd = "git log --pretty=format:'%h - %an, %ar : %s' --since=1.day"
    os.system(cmd+'>'+resultsDir+'gitLog')
	
    os.chdir(cwd)	

#-------------------------------------------------------------------------------
# ModelE specific setup
def setupModelEenv(config, compconfig):
    userconfig =util. ConfigSectionMap(config, 'USERCONFIG')
    branch =  userconfig['repobranch']
    if not branch:
        branch = 'detached'
    resultsDir = userconfig['scratchdir'] + '/results/' + branch
    scratchDir = userconfig['scratchdir'] + '/scratch/' + branch
    makesystem =  userconfig['makesystem']

# the following directories are modelE specific:
    util.mkdir_p(scratchDir+'/decks_repository')
    util.mkdir_p(scratchDir+'/cmrun')
    util.mkdir_p(scratchDir+'/exec')
    util.mkdir_p(scratchDir+'/savedisk')

# We need to get a list of compilers...
    compilers = util.getCompilers(compconfig)
# and libraries/modules info...
    libsconfig = util.ConfigSectionMap(compconfig, 'COMPCONFIG')

# ... to create modelErc file(s) for each compiler
    for comp in compilers:
       if not os.path.exists(scratchDir + comp):
          util.mkdir_p(resultsDir + '/' + comp)
          util.mkdir_p(scratchDir + '/' + comp)
       writeModelErc(libsconfig, scratchDir, comp, makesystem)
       # For planet branches we need SOCRATESPATH in modelErc
       if 'planet' in branch:
	  with open(scratchDir + '/' + comp + '/modelErc.' + comp, 'a') as f:
              f.write('   SOCRATESPATH=/home/damundse/Software/socrates/modele_branch \n')

# Make sure basedirs exist
       if (userconfig['updatebase'] == 'yes'):
		   util.mkdir_p(userconfig['basedir']+'/'+branch+'/'+comp)

#-------------------------------------------------------------------------------
# Write a compiler-specific modelErc file
def writeModelErc(cfg, scratchDir, compiler, makesystem):

   # boos = BUILD_OUT_OF_SOURCE
   boos = 'NO'
   if makesystem == 'makeNew':
      boos = 'YES'

   s = string.Template('\
   DECKS_REPOSITORY=$scr/decks_repository\n\
   CMRUNDIR=$scr/cmrun\n\
   EXECDIR=$scr/exec\n\
   SAVEDISK=$scr/savedisk\n\
   GCMSEARCHPATH=$datadir\n\
   BUILD_OUT_OF_SOURCE=$boos\n\
   COMPILER=$cm\n\
   MPIDISTR=$mn\n\
   MPIDIR=$md\n\
   NETCDFHOME=$nd\n\
   PNETCDFHOME=$pd\n\
   BASELIBDIR5=$bd\n\
   PFUNITSERIALDIR=$p1\n\
   PFUNITMPIDIR=$pn\n\
   OVERWRITE=YES\n\
   OUTPUT_TO_FILES=NO\n\
   VERBOSE_OUTPUT=YES')

   if compiler == 'gfortran':
      modelErc = s.substitute(cm=compiler,\
                              scr=scratchDir,\
                              datadir=cfg['modeldatadir'],\
                              boos=boos,\
                              mn=cfg['gccmpi'],\
                              md=cfg['gccmpidir'],\
                              nd=cfg['gccnetcdf'],\
                              pd=cfg['gccpnetcdf'],\
                              p1=cfg['gccserialpfunitdir'],\
                              pn=cfg['gccmpipfunitdir'],\
                              bd=cfg['gccesmf'])
   elif compiler == 'intel':
      modelErc = s.substitute(cm=compiler,\
                              scr=scratchDir,\
                              datadir=cfg['modeldatadir'],\
                              boos=boos,\
                              mn=cfg['intelmpi'],\
                              md=cfg['intelmpidir'],\
                              nd=cfg['intelnetcdf'],\
                              pd=cfg['intelpnetcdf'],\
                              p1=cfg['intelserialpfunitdir'],\
                              pn=cfg['intelmpipfunitdir'],\
                              bd=cfg['intelesmf'])
   elif compiler == 'nag':
      modelErc = s.substitute(cm=compiler,\
                              scr=scratchDir,\
                              datadir=cfg['modeldatadir'],\
                              boos=boos,\
                              mn=cfg['nagmpi'],\
                              md=cfg['nagmpidir'],\
                              nd=cfg['nagnetcdf'],\
                              pd=cfg['nagpnetcdf'],\
                              p1=cfg['nagserialpfunitdir'],\
                              pn=cfg['nagmpipfunitdir'],\
                              bd=cfg['nagesmf'])
   else:
      modelErc = s.substitute(cm=compiler,\
                              scr=scratchDir,\
                              datadir=cfg['modeldatadir'],\
                              boos=boos,\
                              mn=cfg['gccmpi'],\
                              md=cfg['gccmpidir'],\
                              nd=cfg['gccnetcdf'],\
                              pd=cfg['gccpnetcdf'],\
                              p1=cfg['gccserialpfunitdir'],\
                              pn=cfg['gccmpipfunitdir'],\
                              bd=cfg['gccesmf'])

   rcfile = open(scratchDir + '/' + compiler + '/modelErc.' + compiler, "w")
   rcfile.write(modelErc)
   rcfile.write('\n')
   rcfile.close()
   logger.debug('Created modelErc file for compiler %s', compiler)

#-------------------------------------------------------------------------------
# For in-source builds we need a clone for each rundeck/compiler/mode combo
def setupCloneTasks(config, compconfig, decklist):
    userconfig =util.ConfigSectionMap(config, 'USERCONFIG')
    compilers = util.getCompilers(compconfig)

    cloneTasks = []
    for deck in decklist:
      # since nonProduction rundeck names can be quite long, extract the
      # nonProduction_ part...
        dName = deck.name
        if re.search('nonProduction', deck.name):
            start = deck.name.find('nonProduction') + 14
            dName = deck.name[start:]
         
        for comp in deck.getOpt('compilers').split(','):
            for mode in deck.getOpt('modes').split(','):
                if comp.strip() in compilers:
                   commandString = util.gitCloneCommand(config, dName,
                                                        comp.strip(), mode.strip())
                   cloneTasks.append(commandString)
                else:
                   logger.error('Compiler '+comp+' is not defined in COMPCONFIG')
            
    for t in cloneTasks:
        logger.debug('CLONE TASK %s', t)
    return cloneTasks

#-------------------------------------------------------------------------------
# For out-of source builds, create directory for each rundeck/compiler/mode combo
def setupRuns(config, compconfig, decklist):
    userconfig =util.ConfigSectionMap(config, 'USERCONFIG')
    compilers = util.getCompilers(compconfig)
    scratch = userconfig['scratchdir']
    repo = userconfig['repository']
    branch =  userconfig['repobranch']
    if not branch:
        branch = 'detached'
    repo = scratch + '/scratch/' + branch + '/' + branch
    os.chdir(scratch + '/scratch/' + branch)

    cwd = os.getcwd()
    for comp in compilers:
        util.mkdir_p(cwd+'/'+comp)
        os.chdir(cwd+'/'+comp)
        for deck in decklist:
            dName = deck.name
            if re.search('nonProduction', deck.name):
                start = deck.name.find('nonProduction') + 14
                dName = deck.name[start:]
            for mode in deck.getOpt('modes').split(','):
                adir = dName +  '.' + mode.strip()
                if not os.path.isdir(adir):
                    util.mkdir_p(adir)
    setupModelEenv(config, compconfig)

#-------------------------------------------------------------------------------
# Return a command to submit/execute a [batch] job
def setupScriptTasks(config, compconfig, decklist):
    logger.info('Prepare and execute tasks...')
    compilers = util.getCompilers(compconfig)

    scriptTasks = []
    for deck in decklist:
        for comp in deck.getOpt('compilers').split(','):
            for mode in deck.getOpt('modes').split(','):
                if comp.strip() in compilers:
                    commandString = \
                        createScriptTask(config, compconfig, deck,
                                         comp.strip(), mode.strip())
                    scriptTasks.append(commandString)
                else:
                    logger.error(comp+' is not defined in COMPCONFIG')

    if len(scriptTasks) > 0:
        for t in scriptTasks:
            logger.debug('SCRIPT TASK %s', t)
    else:
            logger.debug('There is nothing to do.')        
    return scriptTasks

#-------------------------------------------------------------------------------
# Creates script to be submitted to batch system OR to be executed interactively
# Batch system is assumed to be the one on NCCS-DISCOVER machines
def createScriptTask(config, compconfig, deck, comp, mode):
    userconfig  = util.ConfigSectionMap(config, 'USERCONFIG')
    modules    = userconfig['modules']
    useBatch   = userconfig['usebatch']
    branch     = userconfig['repobranch']
    if not branch:
        branch = 'detached'
    scriptsDir = userconfig['scriptsdir'] + '/'
    useMods    = userconfig['modules']
    resultsDir = userconfig['scratchdir'] + '/results/' + \
                 branch + '/' + comp
    scratchDir = userconfig['scratchdir'] + '/scratch/' + \
                 branch + '/' + comp
    sponsorID  = userconfig['sponsorid']

    deckName = deck.name
    jobName = deckName
    if re.search('nonProduction', deckName):
       start = deckName.find('nonProduction') + 14
       jobName = deckName[start:]
    filename = resultsDir + '/' + jobName + '.' + mode + '.bash'
    fileHandle = open ( filename, 'w' ) 

    if useBatch == 'yes':
        
        cores = max(deck.getOpt('npes'))

        # If we are just compiling this rundeck
        if deck.getOpt('verification') == 'compileOnly':
            cores = 4
            walltime = '00:10:00'

        # customRun is a 2-month run
        elif deck.getOpt('verification') == 'customRun':
            if re.search('campi', deckName):
                walltime = '4:00:00'
            elif re.search('Tmatrix', deckName):
                walltime = '4:00:00'
            elif re.search('Ttomas', deckName):
                walltime = '4:00:00'
            elif re.search('cadi', deckName):
                walltime = '2:00:00'
            elif re.search('Toma', deckName):
                walltime = '2:00:00'
            elif re.search('obio', deckName):
                walltime = '1:00:00'
            elif re.search('C12', deckName):
                walltime = '0:30:00'
            elif re.search('M20', deckName):
                walltime = '0:30:00'
            else:
                walltime = '1:00:00'
            
        # regular runs (1hr and/or restart)
        else:               
            if 'mpi' in mode: 
                walltime = '01:00:00'

            else: # serial
                cores = 1
                walltime = '00:30:00'
                if re.search('obio', deckName):
                    walltime = '01:00:00'
                elif re.search('cadi', deckName):
                    walltime = '04:00:00'
                elif re.search('Toma', deckName):
                    walltime = '04:00:00'
                elif re.search('lerner', deckName):
                    walltime = '01:00:00'

            # Adjust the walltime for some rundecks
            if re.search('C12', deckName):
                walltime = '00:30:00'
            elif re.search('Mars', deckName):
                walltime = '00:30:00'
            elif re.search('SGP', deckName):
                walltime = '00:10:00'
            elif re.search('campi', deckName):
                walltime = '02:00:00'
            elif re.search('ctomas', deckName):
                walltime = '02:00:00'
            elif re.search('Ttomas', deckName):
                walltime = '03:00:00'
            elif re.search('Tmatrix', deckName):
                walltime = '02:00:00'
            elif re.search('vsd', deckName):
                if 'mpi' in mode: 
                  walltime = '02:00:00' # parallel
                else: 
                  walltime = '03:30:00' # serial

        outname = resultsDir + '/' + jobName + '.' + mode + '.out'
        errname = resultsDir + '/' + jobName + '.' + mode + '.err'
        fileHandle.write ('#!/bin/bash' + '\n')
        fileHandle.write ('#SBATCH -J ' + jobName + '\n')
        fileHandle.write ('#SBATCH -o ' + outname + '\n')
        fileHandle.write ('#SBATCH -e ' + errname + '\n')
        fileHandle.write ('#SBATCH --account='  + sponsorID + '\n')
        fileHandle.write ('#SBATCH --time='     + walltime + '\n')
        fileHandle.write ('#SBATCH --ntasks=' + str(cores) + '\n')
        # Use Haswell NODES
        #fileHandle.write ('#SBATCH --constraint=hasw' + '\n')
        #if walltime == '00:30:00':
        #    fileHandle.write ('#SBATCH --qos=debug' + '\n')
          
    # Create rest of script used in batch OR interactive jobs:

    if modules == 'yes':
        machine = sp.check_output(['uname','-n'])
        # If on NCCS-DISCOVER
        if 'borg' in machine or 'discover' in machine or 'dali' in machine:
            fileHandle.write ('. /usr/share/modules/init/bash' + '\n')
            fileHandle.write ('module purge' + '\n')
            # Need python 2.7.x
            fileHandle.write ('module load other/SSSO_Ana-PyD/SApd_2.1.0' + '\n')
        else:
            fileHandle.write ('#!/bin/bash' + '\n')
            # CC: This is not portable...just my MAC so far
            fileHandle.write ('. /opt/local/share/Modules/3.2.10/init/bash' + '\n')
            fileHandle.write ('module purge' + '\n')

        # Needed due to different naming convention for module names
        compvendor = comp
        if comp == 'gfortran':
            compvendor = 'gcc'

        modsconfig = util.ConfigSectionMap(compconfig, 'COMPCONFIG')
        for mod in modsconfig['modulelist'].split(','):
            if re.search(compvendor, mod):
                for mm in modsconfig[mod].split(','):
                    cmd = 'module load ' + mm +'\n'
                    fileHandle.write (cmd)

    fileHandle.write ('umask 022' + '\n')
    makesystem =  userconfig['makesystem']
    if makesystem == 'makeOld':
        decksDir = scratchDir + '/' + jobName +  '.' + mode + '/decks/'
    else:
        decksDir = scratchDir + '/' + jobName + '.' + mode

    fileHandle.write ('export DECKSDIR=' + decksDir + '\n')
    modelErc = scratchDir + '/modelErc.' + comp
    fileHandle.write ('export MODELERC=' + modelErc + '\n')

    if deck.getOpt('unittest') == 'yes':
        if 'mpi' in mode: 
            cmd = "cat "+modelErc+"| grep PFUNITMPIDIR"+"| awk -F= '{print $2}'"
        else:
            cmd = "cat "+modelErc+"| grep PFUNITSERIALDIR"+"| awk -F= '{print $2}'"
        out = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True)
        pfunitDir = out.communicate()[0].rstrip()
        fileHandle.write ('export PFUNIT=' + pfunitDir + '\n')
        
    fileHandle.write ('cd ' + decksDir + '\n')
    fileHandle.write ('python ' + scriptsDir + '/' + 'regression.py ' + deckName + '\n')
    fileHandle.write (' ' + '\n')
    fileHandle.close()
    
    createRegConfig(config, deck, modelErc, comp, jobName, mode)

    if useBatch == 'yes':
        commandString = 'sbatch ' + filename
    else:
        commandString = 'chmod +x ' + filename + ';' + filename
               
    return commandString

#-------------------------------------------------------------------------------
# Create a config file for regression.py script. 
# Note: there is one config file for each rundeck/compiler combination
def createRegConfig(config, deck, modelErc, comp, jobName, mode):
    cfg  = util.ConfigSectionMap(config, 'USERCONFIG')
    branch     = cfg['repobranch']
    if not branch:
        branch = 'detached'
    resultsDir = cfg['scratchdir'] + '/results/' + \
            branch + '/' + comp
    scratch = cfg['scratchdir'] + '/scratch/' + \
            branch + '/' + comp
    makesystem = cfg['makesystem']
    if makesystem == 'makeOld':
        decksDir = scratch + '/' + jobName + '.' + mode + '/decks/'
    else:
        decksDir = scratch + '/' + jobName + '.' + mode

# Regression tests do not run the regression script in standalone mode    
    standalone = 'no'
    regconfig  = ConfigParser.RawConfigParser()
    regconfig.add_section('regSettings')
    regconfig.set('regSettings', 'rundeck', deck.name)
    regconfig.set('regSettings', 'modelerc', modelErc)
    regconfig.set('regSettings', 'compiler', comp)
    regconfig.set('regSettings', 'modes', mode)
    regconfig.set('regSettings', 'standalone', standalone)
    regconfig.set('regSettings', 'unittest', deck.getOpt('unittest'))
    regconfig.set('regSettings', 'verification', deck.getOpt('verification'))
    regconfig.set('regSettings', 'endtime', deck.getOpt('endtime'))
    npestr = ' '.join(str(e) for e in deck.getOpt('npes'))
    regconfig.set('regSettings', 'npes', npestr)
    regconfig.set('regSettings', 'buildtype', cfg['buildtype'])

    if makesystem == 'makeOld':
        regconfig.set('regSettings', 'repository', cfg['repository'])
    else:
        # Out of source build still pollutes the repository a little bit,
        # specially for nonProduction builds. So, let's make sure we pollute
        # a clone.
        newrepo = cfg['scratchdir']+'/scratch/'+branch+'/'+branch      
        regconfig.set('regSettings', 'repository', newrepo)

    regconfig.set('regSettings', 'branch', branch)
    regconfig.set('regSettings', 'basedir', cfg['basedir'])
    regconfig.set('regSettings', 'updatebase', cfg['updatebase'])
    regconfig.set('regSettings', 'makesystem', makesystem)
    regconfig.set('regSettings', 'resultsdir', resultsDir)
    regconfig.set('regSettings', 'scratchdir', scratch)
    regconfig.set('regSettings', 'decksdir', decksDir)

    filename = decksDir + '/' + deck.name + '.cfg'
    with open(filename, 'w') as configfile:
        regconfig.write(configfile)

#-------------------------------------------------------------------------------
# Create a diff report and notify via email
def sendDiffreport(config, compconfig, eTime):
    userconfig  = util.ConfigSectionMap(config, 'USERCONFIG')
    mailto     = userconfig['mailto']
    branch     = userconfig['repobranch']
    if not branch:
        branch = 'detached'
    resultsDir = userconfig['scratchdir'] + '/results/' + branch
    buildtype  = userconfig['buildtype']
    message    = userconfig['message']
    html       = userconfig['html']
    sortdiff   = userconfig['sortdiff']
    compilers  = util.getCompilers(compconfig)

    diffFile = resultsDir + '/' + 'diffreport.txt'
    fp = open(diffFile, 'w')
    if html == 'yes':
        fp.write('<html><pre>\n')
    fp.write(message + ' \n')
    fp.write('Repository: ' + userconfig['repository'] +  '\n')
    fp.write('-'*80+'\n')
    fp.write('Branch: ' + branch)
    fp.write('  --  Build type: ' + buildtype +  '\n')
    fp.write('-'*80+'\n')
    fp.write('%78s\n' % ('    -REPRODUCIBILITY   '))
    fp.write('%20s%10s%8s%8s%8s%8s%8s%8s\n' % \
        ('RUNDECK', 'COMPILER', 'MODE', 'RUN', 'UNT', 'BAS', 'RST', 'NPE'))
    fp.write('-'*80+'\n')

    if sortdiff == 'yes':
        sp.call('find '+resultsDir+' -name \*.diff -exec cat {} \; >' \
                    +resultsDir + '/' + 'alldiffs', shell=True)
        sp.call('cat '+resultsDir + '/' + 'alldiffs | sort -k 1,1 >' \
                    +resultsDir + '/' + 'sorteddiffs', shell=True)
        with open(resultsDir + '/' + 'sorteddiffs','r') as inf:
            fp.write(inf.read())
        os.remove(resultsDir + '/' + 'alldiffs')
        os.remove(resultsDir + '/' + 'sorteddiffs')
    else:
        for comp in compilers:
            diffs = glob.glob(resultsDir + '/' + comp + '/*.diff')
            for f in diffs:
                with open(f,'r') as inf:
                    fp.write(inf.read())
 
    # In some cases, if a task terminated unexpectedly then the results will
    # not be recorded to a diff file. However, the system should generate an
    # error file (*.err). If such a file exists then we update the results
    # and notify a system error.
    results = [' '*20, ' '*10, ' '*8, 
                        '  -  ', '  -  ', '  -  ', '  -  ', '  -  ']
    for comp in compilers:
        errs = glob.glob(resultsDir + '/' + comp + '/*.err')
        for f in errs:
            if os.path.getsize(f) > 0:
                fname = f.split('/')[-1]
                results[0] = re.split(r'\.(?!\d)', fname)[0]
                results[1] = comp
                results[2] = re.split(r'\.(?!\d)', fname)[1]
                results[3] = 'U'
                util.writeDiff(results, fp)
    
    fp.write('-'*80+'\n')
    hhmmss = time.strftime('%H:%M:%S', time.gmtime(eTime))
    fp.write('Time taken = %s \n' %(hhmmss))
    fp.write('-'*80+'\n')
    fp.write('Legend:\n')
    fp.write('-'*7+'\n')
    fp.write('+   : success\n')
    fp.write('C   : created baseline\n')
    fp.write('Fb  : build failure\n')
    fp.write('F1  : 1hr run-time failure\n')
    fp.write('Fr  : restart run-time failure\n')
    fp.write('F*  : expected failure\n')
    fp.write('U   : unexpected system failure\n')
    fp.write('NUM : number of reproducibility differences\n')
    fp.write('-   : not available\n')
    fp.write('Notes:\n')
    fp.write('-'*6+'\n')
    compconfig = util.ConfigSectionMap(compconfig, 'COMPCONFIG')
    compVers =  compconfig['compiler_versions'].split(",")
    i=0
    for comp in compilers:
        fp.write(comp+' compiler version: '+compVers[i]+'\n')
        i+=1
    fp.write('Results in: ' + resultsDir +  '\n')
    fp.write('-'*80+'\n')
    fp.write( 'Commits from last day:\n')
    with open(resultsDir + '/gitLog', 'r') as inf:
        fp.write(inf.read())
    fp.write( '\n')
    fp.write('-'*80+'\n')
    if html == 'yes':
        fp.write('</pre><html>\n')
    fp.close()

    subject = '"[modelE-regression]" '
    cmd = '/usr/bin/mail -s ' + subject + mailto + ' < ' + diffFile
    if html == 'yes':
        pref = 'mutt -e "set content_type=text/html" -s '
        cmd = pref + subject + mailto + ' < ' + diffFile
    sp.call(cmd, shell=True)

