#Programmer: Nicolas Dullius Mallmann
#Name: control
#Status: Development

import os, re, sys, glob
import libConfig as lc
import threading as th

import inspect

if sys.version_info[0] == 2:
    from Queue import Queue as q
else:
    from queue import Queue as q


def run(configuration_file):

    #Read the configuration file
    conf = lc.Config(configuration_file)

    #Load the list of jobs
    job_order = conf["master"]["jobs_list"]

    #Import all the functions for the jobs
    functions_list = {}
    for job in job_order:
        functions_list[job] = __import__(name = conf[job]["module"], globals = globals(), locals = locals()).__dict__[conf[job]["function"]]

    #Test the existence of the root_path option in the configuration file
    if not conf.has_option("master", "root_path"):
        print("Warning! 'root_path' was not set. Setting it to './Results'!")
        conf.set("master", "root_path", "./Results")

    #Create the root_path directory if it doesn't exist
    root_path = conf["master"]["root_path"]
    if not os.path.exists(root_path):
        os.makedirs(root_path)

    #Get the initial data
    initial_is_list = conf["master"]["list"]
    data_ref = conf["master"]["initial_data"]

    #If the file is a list of files
    if initial_is_list:
        initial_data_list = []
        for data in open(data_ref).readlines():
            data = re.split("\s*#", data.rstrip("\n"))[0]
            if len(data) > 0:
                initial_data_list.append(data)

    #If not a list, test if the initial data is a directory or a single file
    else:
        #If initial data is a directory, search for all files inside it (with the chosen extension)
        if os.path.isdir(data_ref):
            extension = conf["master"]["initial_data_extension"].lstrip(".")
            initial_data_list = glob.glob(os.path.join(data_ref, "*{0}".format(extension)))
        #If initial data is a file
        elif os.path.isfile(data_ref):
            initial_data_list = [data_ref]
        else:
            raise RuntimeError("Initial data is not valid.")

    #Run n procresses at the same time
    number_of_processes = conf["master"]["processes"]

    if not isinstance(number_of_processes, int):
        print("Warning! 'processes' (in the master category) must be an integer parameter! Setting it to 1!")
        conf.set("master", "processes", 1)
        number_of_processes = 1
    elif number_of_processes <= 0:
        print("Warning! 'processes' (in the master category) must be an integer greater than zero! Setting it to 1!")
        conf.set("master", "processes", 1)
        number_of_processes = 1
    elif number_of_processes > len(initial_data_list):
        print("Warning! 'processes' (in the master category) is higher than the number of initial data! Setting it to {0}!".format(len(initial_data_list)))
        conf.set("master", "processes", len(initial_data_list))
        number_of_processes = len(initial_data_list)

    #Initialize queue of processes
    list_of_processes = []
    process_queue = q(number_of_processes)

    #Start the distribuition of jobs
    for data in initial_data_list:
        process_queue.put(1)

        new_config = conf.copy()
        new_config.set(job_order[0], "input", data)

        list_of_processes.append(process(job_order, functions_list, new_config, process_queue))

        list_of_processes[-1].start()

        #Wait for the rest of the processes to end
        process_queue.join()

class process(th.Thread):
    def __init__(self, job_order, functions_list, configuration_object, queue):
        th.Thread.__init__(self)
        self.job_order = job_order
        self.functions_list = functions_list
        self.conf = configuration_object
        self.queue = queue

    def run(self):
        self.queue.get()

        n = 0
        for job in self.job_order:
            # Run the function passing the necessary parameters
            function = self.functions_list[job]
            function(self.conf, n, self.job_order)
            n = n+1

        #Tell the program that the job chain is done
        self.queue.task_done()

if __name__ == "__main__":

    if len(sys.argv) == 2:
        run(sys.argv[1])
    elif len(sys.argv) == 1:
        print("Missing configuration file!\n")
    else:
        print("Too many arguments, using just the first!")
        run(sys.argv[1])
    
