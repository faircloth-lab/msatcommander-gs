import time
import multiprocessing
import subprocess


def q_runner(n_procs, list_item, function, *args):
    '''generic function used to start worker processes'''
    task_queue      = multiprocessing.Queue()
    results_queue   = multiprocessing.JoinableQueue()
    if args:
        arguments = (task_queue, results_queue,) + args
    else:
        arguments = (task_queue, results_queue,)
    results = []
    # reduce processer count if proc count > files
    if len(list_item) < n_procs:
        n_procs = len(list_item)
    for l in list_item:
        task_queue.put(l)
    for _ in range(n_procs):
        p = multiprocessing.Process(target = function, args = arguments).start()
    for _ in range(len(list_item)):
        # join the queue until we're finished processing
        results.append(results_queue.get()) 
        print results
        results_queue.task_done()
    #tell child processes to stop
    for _ in range(n_procs):
        task_queue.put('STOP')
    results_queue.join()
    return results

def worker1(input, output):
    for c in iter(input.get, 'STOP'):
        # do some stuff
        two_bit = subprocess.Popen('sleep 2', shell=True, stdout=subprocess.PIPE, stderr = subprocess.PIPE).communicate(None)
        # stick the ouput somewhere
        output.put(str(c) + 'a')



n_procs = 1
list_item = range(10)
function  = worker1
q_runner(n_procs, list_item, function)