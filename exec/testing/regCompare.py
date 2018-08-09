import os
import subprocess as sp
import logging
import regUtils as utils

"""
  Compare model results with those in the baseline location.
  If no baseline location is specified then comparison will be skipped.
"""
def base(run, endTime, npes=1):
    logger = logging.getLogger('COMPBAS ')
    logger.info('Compare '+str(endTime)+'hr and base run')

    # Skip comparison if no baseDir location is given
    if run.baseDir == '.':
        logger.info('No baseline directory - nothing to do')
        return

    if run.makesystem == 'makeOld':
        prefix = run.name + '/'
    else:
        prefix = run.decksDir+'/'+run.name+'/'
        
    # Check if checkpoint file exists:
    current = prefix + \
        utils.checkpointName(run.name, run.mode, str(endTime)+'hr', npes)
    if not os.path.exists(current):
        logger.error('CHECKPOINT file '+current+' does not exist.')
        return

    # Check if baseline file exists:
    theBase = run.baseDir + '/' + run.shortName + '.' + str(endTime)+'hr'
    tmpBase = theBase
    # If baseline does not exist then copy it to baseline directory
    if not os.path.exists(theBase):
        logger.warning('Baseline file '+theBase+' does not exist!')
        if run.updateBase == 'yes':
            if utils.copyFile(current, theBase) == 0:
                logger.info('Created BASELINE')
                run.results[5] = run.createMark
            else:
                logger.error('Error in: cp '+current+' ' +theBase)
        else:
            run.results[5] = run.naMark

    # If baseline exists then execute diffreport
    else:
        logger.debug('Compare '+current+' '+theBase)
        n = getNumDiffs(run, theBase, current)
        # No differences
        if n == '0':
            run.results[5] = run.successMark
        # If there are differences then update the baseline dir (if specified)
        # and keep a (temporary) copy of the old one
        else:
            run.results[5] = '{: ^5}'.format(n)
            logger.warning('Baseline reproducibility failed!!!')
            if run.updateBase == 'yes':
                if utils.copyFile(current, theBase) == 0:
                    logger.info('Replaced BASELINE')
                else:
                    logger.error('Error in: cp '+current+' ' +theBase)
                # TODO: Make this configurable
                #timestr = time.strftime("%Y%m%d-%H%M%S")
                #if utils.copyFile(theBase, tmpBase+'-'+timestr) == 0:
                #    logger.info('Saved old BASELINE')
                #else:
                #    logger.error('Error in: cp '+theBase+' ' +tmpBase+'-'+timestr)
            else:
                logger.info('BASELINE was not replaced: updateBase=' \
                            + run.updateBase)

"""
  Compare continuous-run vs restart run
"""
def restart(run, npes=1):
    logger = logging.getLogger('COMPRST ')
    logger.info('Compare '+str(run.endTime)+'hr and restart run')
    prefix = run.name + '/'

    # Check if continuous run result exists:
    fileCON = prefix + \
        utils.checkpointName(run.name, run.mode, str(run.endTime)+'hr', npes)
    if not os.path.exists(fileCON):
        logger.error('CHECKPOINT file '+fileCON+' does not exist.')
        return

    # Check if restart run result exists:
    fileRST = prefix + \
        utils.checkpointName(run.name, run.mode, 'restart', npes)
    if not os.path.exists(fileRST):
        logger.error('RESTART file '+fileRST+' does not exist.')
        return

    logger.debug('Compare '+fileCON+' '+fileRST)
    n = getNumDiffs(run, fileCON, fileRST)
    if n == '0':
        run.results[6] = run.successMark
    else:
        # SCM rundeck is not restart reproducible
        if 'SGP' in run.name:
            run.results[6] = run.failMark+'*'
        else:
            run.results[6] = '{: ^5}'.format(n)
            logger.warning('Restart reproducibility failed!!!')

"""
  Compare SERIAL vs MPI
"""
def nPE(runMPI, endTime, npes):
    logger = logging.getLogger('COMPNPE ')
    logger.info('Compare serial and '+runMPI.mode+ \
                    ' (' +str(npes)+' npes) runs')

    if runMPI.makesystem == 'makeOld':
        ddir = '/decks/'
    else:
        ddir = '/'

    # Check if SERIAL result exists:
    fileSER = runMPI.savedisk + '/' + \
        runMPI.shortName+'.serial.' + \
        runMPI.compiler + '/' + \
        runMPI.shortName+'.serial.' + \
        runMPI.compiler + '.' + str(endTime) + 'hr'
    if not os.path.exists(fileSER):
        logger.warning('SERIAL file '+fileSER+' does not exist!')
        return

    # Check if MPI result exists:
    fileMPI = runMPI.name + '/' + \
       utils.checkpointName(runMPI.name, runMPI.mode, str(endTime)+'hr', npes)
    if not os.path.exists(fileMPI):
        logger.warning('MPI file '+fileMPI+' does not exist!')
        return

    logger.debug('Compare '+fileSER+' '+fileMPI)
    n = getNumDiffs(runMPI, fileSER, fileMPI)
    if n == '0':
        runMPI.results[7] = runMPI.successMark
    else:
        runMPI.results[7] = '{: ^5}'.format(n)
        logger.warning('NPE reproducibility failed!!!')
  
"""
  Get number of -diffs- when running diffreport
"""
def getNumDiffs(run, file1, file2):
    # Create a file to hold doffreport output
    file1str = file1.split('/')
    file2str = file2.split('/')
    diffFile = run.resultsDir+'/'+file1str[-1]+'-vs-'+file2str[-1]
    # Command to run diffreport executable
    cmd1 = getDiffexe()+' '+file1+' '+file2+' > '+diffFile
    diff1 = sp.Popen(cmd1, shell=True)
    diff1.wait()
    # Command to count number of diffs in diffFile
    cmd2 = 'cat '+diffFile+' | grep diffs | wc -l'
    diff2 = sp.Popen(cmd2, stdout=sp.PIPE, shell=True)
    numDiffs = diff2.communicate()[0]
    diff2.wait()
    # If there are no diffreport differences remove the diffFile
    if os.stat(diffFile).st_size == 0:
        os.remove(diffFile)
    return ''.join(numDiffs.split())

"""
   Get diffreport executable path
"""
def getDiffexe():
    # This is needed to find diffreport.x, assumed to be in $HOME/bin
    os.environ["PATH"] += os.pathsep + os.environ["HOME"] + '/bin'
    diffreportExe = utils.which('diffreport.x')
    if diffreportExe is None:
        print 'No available diffreport.x. Will use diff'
        diffreportExe = 'diff'
    return diffreportExe
