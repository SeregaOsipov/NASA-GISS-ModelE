# This module contains functions that help setup, build and run modelE rundecks
import sys
import os
import re
import subprocess as sp
import logging
import ConfigParser
import regUtils as utils
import regTest as t

""" 
  This class assigns settings used to perform regression tests on a rundeck 
"""
class newRundeck(t.regTest):

    def __init__(self, sourceName='nonProduction_E_AR5_C12'):
        # Overide name
        if sourceName:
            self.name      = sourceName
        else:
            self.name      = 'nonProduction_E_AR5_C12'
        super(newRundeck, self).__init__(sourceName)

        # modelE rundeck defaults    
        self.standalone    = 'yes'
        self.branch        = ''
        self.updateBase    = 'no'
        self.baseDir       = '.'
        self.repository    = '..'
        self.decksDir      = '.'
        self.makesystem    = 'makeOld'
        self.compiler      = '' 
        self.savedisk      = '.'
        # Now override defaults with values specified in config file
        getConfiguration(self)

    #  Setup a logging object for each rundeck. Note that the output file gets
    #  logging DEBUG while STDOUT only gets logging INFO in order to minimize 
    #  verbosity.
    def setLogging(self):
        # Note filemode is 'append' because we want one logger per rundeck
        logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename=self.resultsDir+'/'+self.name
                             +'-regression.log',
                    filemode='a')
        stdoutLog = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter('%(name)s : %(message)s')
        stdoutLog.setFormatter(formatter)
        stdoutLog.setLevel(logging.INFO)
        logger = logging.getLogger()
        logger.addHandler(stdoutLog)
        

"""
  This class defines, among other things, a rundecks's run time options
  passed to the the makefile and its build/run methods.
  Note that a "run" is a "rundeck".
"""
class newRun(newRundeck):

    def __init__(self, rundeck, mode):
        super(newRun, self).__init__(rundeck.name)
        self.runsrc = rundeck.name
        # shortName is introduced to remove the nonProduction part of a
        # rundeck name
        shortName = rundeck.name
        if 'nonProduction' in rundeck.name:
            shortName = rundeck.name[14:]
        # A RUN name is composed of rundeck name, mode and compiler
        self.shortName = shortName
        self.name = shortName+'.'+mode+'.'+rundeck.compiler
        # The following are (cmd)strings used by gnu-make
        self.runCmd = 'RUN='+self.name
        self.runSrcCmd = 'RUNSRC='+rundeck.name
        if mode == 'serial':
            self.mode = 'serial'
            self.modeCmd = 'MPI=NO'
        else:
            self.mode = 'mpi'
            self.modeCmd = 'MPI=YES'
        self.xflags = ' '
        if rundeck.buildType == 'debug':
            flags='"-O0 -g"'
            self.xflags += 'EXTRA_FFLAGS='+flags
        elif rundeck.buildType == 'traps':
            self.xflags += 'COMPILE_WITH_TRAPS=YES'
        # Run results are stored in an array
        self.results = [shortName, rundeck.compiler, mode, 
                        '  -  ', '  -  ', '  -  ', '  -  ', '  -  ']
        self.successMark   = '+'
        self.failMark      = 'F'
        self.createMark    = 'C'
        self.naMark        = '-'
        # In case we want to debug this run:
        if os.environ.has_key('DEBUG'):
            self.debug = os.environ['DEBUG']
        else:
            self.debug = False
        # modelE return code for a successful run:
        self.OKrc = '13'

    # ------------------------------------------------------
    #  Make system calls and record results of a regression test
    def sysCmd(self, cmdStr, resIdx, stageID):

        logger = logging.getLogger('SYSTEM  ')
        status = 1
        numfail = 0
        if self.debug:
            logger.info(cmdStr)
        else:
            logger.debug(cmdStr)
            if stageID == 't':
                if os.environ.has_key('PFUNIT'):
                    # --- Run make tests command ---
                    proc = sp.Popen(cmdStr, stdout=sp.PIPE, stderr=sp.PIPE,
                                    shell=True)
                    # --- Check output "interactively"
                    while True:
                       output = proc.stdout.readline()
                       if output == '' and proc.poll is not None:
                           break
                       if output:
                           logger.debug(output.strip())
                           if 'Failures' in output:
                               outlist =  output.split(' ')
                               # Extract number of unit test failures, if any
                               numfail = int(re.search(r'\d+', outlist[4]).group())                            
                    proc.wait()
                    # --- Done with tests ---
                    status = proc.returncode
                else:
                    logger.error('PFUNIT environment variable has not been set.')
            else:
                makeLog = self.resultsDir + '/'  + self.name + '-make.log'
                with open(makeLog,'a') as f:
                    proc = sp.Popen(cmdStr,stdout=f,stderr=f,shell=True)
                    proc.wait()
                    status = proc.returncode
                                 
            logger.debug('Return code: ' + str(status))

            if (status == 0):
                self.results[resIdx] = self.successMark
            else:
                if stageID == 't':
                    self.results[resIdx] = str(numfail)
                else:
                    logger.error(cmdStr+': FAILED')
                    self.results[resIdx] = self.failMark+stageID
                raise RuntimeError(cmdStr+': failed')
    
    # ------------------------------------------------------
    def build(self):
        logger = logging.getLogger('BUILD   ')
        logger.info(self.name + ' ' + self.modeCmd + ' ' + self.xflags)

        stageID = 'b'

        if self.makesystem == 'makeOld':
            try:
                cmd = 'make rundeck '+self.runCmd+' '+self.runSrcCmd
                self.sysCmd(cmd, 3, stageID)
            except RuntimeError, e:
                return 1

            if self.standalone == 'yes':
                try:
                    cmd = 'make --quiet clean'
                    self.sysCmd(cmd, 3, stageID)
                except RuntimeError, e:
                    return 1
            try:
                cmd = 'make -j gcm '+self.runCmd+' '+self.modeCmd+' '+self.xflags
                self.sysCmd(cmd, 3, stageID)
            except RuntimeError, e:
                return 1
        else:
            try:
                configure = self.repository+'/exec/configure '
                cmd = configure+self.name+' '+self.runsrc
                self.sysCmd(cmd, 3, stageID)
            except RuntimeError, e:
                return 1
            try:
                cmd = 'make -j gcm'
                self.sysCmd(cmd, 3, stageID)
            except RuntimeError, e:
                return 1

        if self.unitTest == 'yes':
            logger.info('Run unit tests...')
            try:
                cmd = 'make tests '+self.runCmd+' '+self.modeCmd
                self.sysCmd(cmd, 4, 't')
            except RuntimeError, e:
                return 1
        return 0
           
    # ------------------------------------------------------
    # Sets up and runs a 1hr simulation
    def hrRun(self, npes=1):
        logger = logging.getLogger('RUN1HR  ')
        logger.info(self.name + ', ' + self.mode + ', npes=' + str(npes))

        stageID = '1'

        if self.makesystem == 'makeOld':
            try:
                cmd =  'make -j setup '+self.runCmd+' '+self.modeCmd+' '+self.xflags
                self.sysCmd(cmd, 3, stageID)
            except RuntimeError, e:
                return 1
        else:
            try:
                self.sysCmd('make -j setup', 3, stageID)
            except RuntimeError, e:
                return 1

        try:
            rune = self.repository+'/exec/runE '
            cmd = rune+self.name+' -np '+str(npes)+' -cold-restart'  
            self.sysCmd(cmd, 3, stageID)
        except RuntimeError, e:
            return 1

        try:
            cmd = 'cd '+self.name+ '; test `head -1 run_status` -eq ' + self.OKrc
            self.sysCmd(cmd, 3, stageID)
        except RuntimeError, e:
            return 1

        try:
            cmd = 'cd ' + self.name + '; cp fort.2.nc ' \
                + utils.checkpointName(self.name, self.mode, '1hr', npes)
            self.sysCmd(cmd, 3, stageID)
        except RuntimeError, e:
            return 1

        logger.info(self.name + ' is DONE')
        return 0

    # ------------------------------------------------------
    # Run up to ENDTIME hrs with checkpoint at ENDTIME-1 hrs
    def restartRun(self, npes=1, endTime=25):
        logger = logging.getLogger('RUNRST  ')

        restart = './'+self.name
        if self.mode == 'mpi':
            restart += ' -np ' + str(npes)

        logger.info(self.name + ', ' + self.mode + ', npes=' + str(npes) + \
                        ', endTime=' + str(endTime))

        stageID = 'r'
        checkPt = endTime - 1
        ndisk = checkPt * 2
        if endTime > 24:
            newTime = ' ' + str(ndisk) + ' 2 1'
        else:
            newTime = ' ' +  str(ndisk) + ' 1 ' + str(endTime)

        try:
            cmd = self.repository+'/exec/editRundeck.sh '+self.name+newTime 
            self.sysCmd(cmd, 3, stageID)
        except RuntimeError, e:
            return 1

        # Run make setup again to update the rundeck-derived I file
        if self.makesystem == 'makeOld':
            try:
                cmd =  'make -j setup '+self.runCmd+' '+self.modeCmd+' '+self.xflags
                self.sysCmd(cmd, 3, stageID)
            except RuntimeError, e:
                return 1
        else:
            try:
                self.sysCmd('make -j setup', 3, stageID)
            except RuntimeError, e:
                return 1

        try:
            rune = self.repository+'/exec/runE '
            cmd = rune+self.name+' -np '+str(npes)+' -cold-restart'  
            self.sysCmd(cmd, 3, stageID)
        except RuntimeError, e:
            return 1

        try:
            cmd = 'cd ' + self.name + '; cp fort.1.nc ' \
                +utils.checkpointName(self.name, self.mode, str(endTime)+'hr', npes)
            self.sysCmd(cmd, 3, stageID)
        except RuntimeError, e:
            return 1

        try:
            cmd = 'cd ' + self.name + '; cp fort.2.nc fort.1.nc; rm -f run_status'
            self.sysCmd(cmd, 3, stageID)
        except RuntimeError, e:
            return 1 

        try:
            cmd = 'cd ' + self.name + '; ' + restart \
                + '; test `head -1 run_status` -eq ' + self.OKrc
            self.sysCmd(cmd, 3, stageID)
        except RuntimeError, e:
            return 1

        try:
            cmd = 'cd ' + self.name + ';cp fort.2.nc ' \
                + utils.checkpointName(self.name, self.mode, 'restart', npes)
            self.sysCmd(cmd, 3, stageID)
        except RuntimeError,e:
            return 1

        # Reset INPUTZ settings (changed by editRundeck.sh) for next MPI run
        if len(self.npes) > 1:
            try:
                makecmd = 'make rundeck '+self.runCmd+' '+self.runSrcCmd
                self.sysCmd(makecmd, 3, stageID)
            except RuntimeError, e:
                return 1

        logger.info(self.name + ' is DONE')
        return 0

   
    # ------------------------------------------------------
    # Run 2 month run
    def customRun(self, npes=1):
        logger = logging.getLogger('RUN2MOS ')
        logger.info(self.name + ', ' + self.mode + ', npes=' + str(npes) \
                    + ', 2 month run')

        stageID = 'L'

        # Most rundecks, except hycom, have MONTHI=12, so MONTHE=2 
        # For hycom, MONTHI=1, so set MONTHE=3. Note ndisk=480.
        if 'h4c' in self.name:
            newTime = ' 480 1 0 3'
        else:
            newTime = ' 480 1 0 2'

        try:
            self.sysCmd(self.repository + '/exec/editRundeck.sh ' + self.name + \
                           newTime, 3, stageID)
        except RuntimeError, e:
            return 1

        if self.makesystem == 'makeOld':
            try:
                self.sysCmd('make -j setup ' + self.runCmd + ' ' + self.modeCmd + \
                               ' ' + self.xflags, 3, stageID)
            except RuntimeError, e:
                return 1
        else:
            try:
                self.sysCmd('make -j setup', 3, stageID)
            except RuntimeError, e:
                return 1

        try:
            self.sysCmd(self.repository + '/exec/runE ' + self.name + ' -np ' + \
                           str(npes) + ' -cold-restart', 3, stageID)
        except RuntimeError, e:
            return 1

        logger.info(self.name + ' is DONE')
        return 0


# ------------------------------------------------------
# Get rundeck's configuration
def getConfiguration(rundeck):

        if os.environ.has_key('MYCONFIGDIR'):
            myConfigDir = os.environ['MYCONFIGDIR']
        else:
            myConfigDir = '.'
        configfile = myConfigDir + '/' + rundeck.name + '.cfg'
    
        if os.path.isfile(configfile):
            config = ConfigParser.ConfigParser()
            config.read(configfile)
            # rundeck settings
            rundeck.name      = config.get('regSettings', 'rundeck')
            rundeck.modelerc  = config.get('regSettings', 'modelerc')
            rundeck.compiler  = config.get('regSettings', 'compiler')
            # rundeck.modes is a list - config setting is a string, so...
            rundeck.modes     = config.get('regSettings', 'modes').split(',')
            # npes is a string in the config file but a list of ints 
            # in the data structure...
            npes = config.get('regSettings', 'npes')
            nint = []
            for n in npes.split():
                nint.append(int(n))
            rundeck.npes = nint
            rundeck.standalone = config.get('regSettings', 'standalone')
            rundeck.verification = config.get('regSettings', 'verification')
            rundeck.unitTest = config.get('regSettings', 'unittest')
            rundeck.endTime = int(config.get('regSettings', 'endtime'))
            rundeck.buildType = config.get('regSettings', 'buildtype')
            # system settings
            rundeck.baseDir = config.get('regSettings', 'basedir')
            rundeck.branch = config.get('regSettings', 'branch')
            rundeck.updateBase = config.get('regSettings', 'updatebase')
            rundeck.resultsDir = config.get('regSettings', 'resultsdir')
            rundeck.scratchDir = config.get('regSettings', 'scratchdir')
            rundeck.decksDir = config.get('regSettings','decksdir')
            rundeck.repository = config.get('regSettings', 'repository')
            rundeck.makesystem = config.get('regSettings', 'makesystem')
        else:
            print ' *** Configuration file does not exist *** ' + configfile
            sys.exit(1)

        # Extract SAVEDISK from modelErc file - needed in compareNPE            
        cmd = "cat "+rundeck.modelerc+"| grep SAVEDISK"+"| awk -F= '{print $2}'"
        out = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True)
        rundeck.savedisk = out.communicate()[0].rstrip()

        for mode in rundeck.modes:
            if mode != 'serial' and mode != 'mpi':
                print ' *** Incorrect mode *** ' + mode

        # This avoid errors in compareBase() when baseDir is not available
        if rundeck.baseDir != '.':
			# For detached branches basedir is associated with master branch
            if rundeck.branch == 'detached':
                rundeck.baseDir =  rundeck.baseDir + '/master/' + rundeck.compiler	
            else:
                rundeck.baseDir =  rundeck.baseDir + '/' + rundeck.branch + '/' \
                + rundeck.compiler

