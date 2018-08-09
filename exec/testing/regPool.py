import time
import subprocess
import shlex
import multiprocessing
import logging

logger = logging.getLogger('pool')

#-------------------------------------------------------------------------------
class Worker(multiprocessing.Process):
    def __init__(self, task_queue, result_queue):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                logger.debug('%s: Exiting' % proc_name)
                self.task_queue.task_done()
                break
            logger.debug('%s: %s' % (proc_name, next_task))
            answer = next_task()
            self.task_queue.task_done()
            self.result_queue.put(answer)
        return

#-------------------------------------------------------------------------------
class Task(object):
    def __init__(self, a):
        self.a = a

    def __call__(self):
        rc =  subprocess.Popen(self.a, shell=True)
        while rc.poll() is None:
            time.sleep(10)
        if rc.returncode !=0:
            logger.debug('%r failed: %s' % (self.a, rc))
        logger.debug('%r is done' % (self.a))

    def __str__(self):
        return '%s' % (self.a)

#-------------------------------------------------------------------------------
class Batch(object):
    jobs = []
    def __init__(self, cmd):
        self.cmd = cmd

    def __call__(self):
        self.slurmCommand()
        while True:
            # If list is empty then we are done
            if len(self.jobs) == 0:
                break                      
            for job in self.jobs:
                rc = subproc('squeue -j '+str(job)+' -t PD,R -h -o %t')
                # If job is done, remove from list
                if not rc:
                    logger.debug( '...' + job + ' is done')
                    self.jobs.remove(job)
                else:
                    if 'PD' in rc:
                        logger.debug( '...' + job + ' is pending')
                    elif 'R' in rc:
                        logger.debug( '...' + job + ' is running')
                    elif 'U' in rc:
                        logger.debug( '...' + job + ' is terminating')
                    else:
                        pass
                    time.sleep(30)

    def __str__(self):
        return '%s' % (self.cmd)

    # This function submits a batch job under SLURM and creates a jobs list
    # NCCS-DISCOVER only.
    def slurmCommand(self):
        rc = subproc(self.cmd % vars())
        logger.debug(rc)
        # output example: "Submitted batch job 12345", so we
        # parse the output from rc and append job ID to jobs list
        self.jobs.append(shlex.split(rc)[3])

#-------------------------------------------------------------------------------
def subproc(cmd):
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    return out

#-------------------------------------------------------------------------------
def runCommands(cmds, useBatch):
    # Establish communication queues
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()
    
    # Start workers
    num_jobs = len(cmds)
    logger.debug( 'Creating %d workers' % num_jobs)
    workers = [ Worker(tasks, results)
                  for i in xrange(num_jobs) ]
    for w in workers:
        w.start()
    
    # Enqueue jobs
    for cmd in cmds:
        if useBatch == 'yes':
            tasks.put(Batch(cmd))
        else:
            tasks.put(Task(cmd))

    # add one stop value per worker to the job queue
    for i in xrange(num_jobs):
        tasks.put(None)

    # Wait for all of the tasks to finish
    tasks.join()
  





  


