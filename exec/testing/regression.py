"""
  This script executes modelE serial and MPI run combinations and generates 
  results used to verify various reproducibility measures.

  The script can be executed from the decks directory as follows:

     ../exec/testing/regression.sh [RUNSRC1 RUNSRC2 ... RUNSRCN -r]

  See the documentation in regression.sh 

  If executed standalone as follows:
 
      python  ../exec/testing/regression.py  <runsource1> [<runsource2> ...]
 
  then this requires a configuration file name <runsource>.cfg for each runsource.
  When called from the higher level driver (reg), this script's results
  are used to generate a report of regression testing reproducibility checks.
  
"""

import sys
import os
import logging
import ConfigParser
import regUtils as util
import regCompare as comp
import regRuns as runsrc

#-------------------------------------------------------------------------------
def compare(runs, fileH):
    for run in runs:
        
        # Internal consistency checks
        if run.mode == 'serial':
            if run.verification == 'restartRun':
                comp.restart(run, run.endTime)
        else:
            for npes in run.npes:
                if run.verification == 'restartRun':
                    comp.restart(run, npes=npes)
                # Compare NPE vs serial
                comp.nPE(run, run.endTime, npes)
            
        # Baseline checks
        if run.mode == 'serial':
            comp.base(run, 1) # 1hr run
            if run.verification == 'restartRun':
                comp.base(run, run.endTime)
        else:
            for npes in run.npes:
                comp.base(run, 1, npes=npes) # 1hr run
                if run.verification == 'restartRun':
                    comp.base(run, run.endTime, npes=npes)
                    
        util.writeDiff(run.results, fileH)
    
#-------------------------------------------------------------------------------
def main():
    OK = 0

    # List of runSources to verify specified on command line
    runSources = []
    if len(sys.argv) > 1:
        for arg in sys.argv[1:]:
            runSources.append(arg)
    else:
        runSources.append('nonProduction_E_AR5_C12')

    for source in runSources:

        src = runsrc.newRundeck(source)

        # Setup a logging stream for this rundeck
        src.setLogging()
        logger = logging.getLogger('MAIN    ')
        
        # Create a list of rundeck run configurations for each mode
        runs = []
        for mode in src.modes:
            runs.append(runsrc.newRun(src, mode))

        # All the work is done from the modelE decks directory
        os.chdir(src.decksDir)

        diffFile = src.resultsDir + '/' + src.name + '.diff'
        fileH = open(diffFile, 'a')

        logger.info('Verifying ' + src.name + ': ' + src.verification)

        for run in runs:

            serBuildResult = OK
            mpiBuildResult = OK

            if run.mode == 'serial':
                serBuildResult = run.build()
                if src.verification == 'compileOnly':
                    continue
                if serBuildResult == OK:
                    # Always run 1hr
                    rc = run.hrRun()
                    if rc == OK and src.verification == 'restartRun':
                        rc = run.restartRun(endTime=run.endTime)
                    else:
                        continue
                else:
                    continue
            else: # MPI
                mpiBuildResult = run.build()
                if run.verification == 'compileOnly':
                    continue
                if mpiBuildResult == OK:
                    if run.verification == 'customRun':
                        for npes in run.npes:
                            rc = run.customRun(npes=npes)
                            if rc != OK:
                                continue
                    else:    
                        for npes in run.npes:
                            # Always run 1hr
                            rc = run.hrRun(npes=npes)
                            if rc == OK and run.verification == 'restartRun':
                                rc = run.restartRun(npes=npes, endTime=run.endTime)
                            else:
                                continue
                else:
                    continue

            logger.info(src.name + ' ' + run.mode + ' runs complete.')
            
        compare(runs, fileH)
        
        fileH.close()

        logger.info(src.name + ' verification complete.')


    logger.info('Regression testing is done.')            
                    
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    main()

                   
