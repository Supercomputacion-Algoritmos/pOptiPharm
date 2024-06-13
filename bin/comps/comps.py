from email import message
from math import ceil, floor
import subprocess
import sys
import pandas as pd
import statistics
import os
import glob
import threading

db_opini = "FMOL-JERO"
db_location = "/home/jero/uni/DB/fda/" # Absolute path. Include last '/'
# db_query = "DB00460.mol2" # Must include extension
db_query = sys.argv[1]
db_tries = 2
db_threads = [1, 2, 4, 8]
server_threads = 96

def sec_run(cmd, opini, query, target):
    sec_run = subprocess.Popen(cmd.format(opini, query, target), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, universal_newlines=True)
    sec_out, sec_err = sec_run.communicate()
    sec_lines = sec_out.splitlines()
    sec_data = sec_lines[0].split()

    return float(sec_data[4]), float(sec_data[16])

def par_run(cmd, opini, query, target, threads, pinned):
    par_run = subprocess.Popen(cmd.format(opini, query, target, threads, pinned), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, universal_newlines=True)
    par_out, par_err = par_run.communicate()
    par_lines = par_out.splitlines()
    par_data = par_lines[len(par_lines) - 2]
    par_data = par_data.split()

    return float(par_data[4]), float(par_data[16])


mutex = threading.Lock() # Table mutex

header = ["OP Version", "Query", "Target", "Pinned", "Threads", "Avg. Time (s)", "Std. Dev. Time (s)", "Avg. Value", "Std. Dev. Value", "Relative speed-up"]
table = pd.DataFrame(columns=header)

sec_cmd = "./OPShapeSimilarity -c {} -q {} -d {} -h 1 -S"
par_cmd = "./pOPShapeSimilarity -c {} -q {} -d {} -h 1 -S -th {} -pn {}"

output_csv = "./{}-vs-all.csv".format(db_query.split(".")[0])

def measure(mine: list):
    global table
    for target in mine:
        db_target = os.path.basename(target)
        print("[Thread {}]Comparing {} vs {}.".format(threading.current_thread().name, db_query, db_target))

        # Skip itself
        if db_target == db_query:
            continue

        tmp_table = pd.DataFrame(columns=header)

        # Secuential OP run
        sec_times = list()
        sec_values = list()
        for t in range(db_tries):
            sec_time, sec_value = sec_run(sec_cmd, db_opini, db_location + db_query, target)
            sec_times.append(sec_time)
            sec_values.append(sec_value)

        tmp_table.loc[len(tmp_table.index)] = ["Secuential", db_query, db_target, "-", 1, statistics.mean(sec_times), statistics.stdev(sec_times), statistics.mean(sec_values), statistics.stdev(sec_values), "-"]
        
        # Threading OP run
        for th in db_threads:
            par_times = list()
            par_values = list()
            
            for t in range(db_tries):
                par_time, par_value = par_run(par_cmd, db_opini, db_location + db_query, target, th, 0)
                par_times.append(par_time)
                par_values.append(par_value)

            tmp_table.loc[len(tmp_table.index)] = ["Threaded", db_query, db_target, 0, th, statistics.mean(par_times), statistics.stdev(par_times), statistics.mean(par_values), statistics.stdev(par_values), statistics.mean(sec_times) / statistics.mean(par_times)]
    
        mutex.acquire()
        table = pd.concat([table, tmp_table], ignore_index=True)
        table.to_csv(output_csv, index=False)
        mutex.release()



db_all_mols = glob.glob(db_location + "*.mol2")
db_all_size = floor(len(db_all_mols) / floor(server_threads / 10))
db_all_rem = len(db_all_mols) % floor(server_threads / 10)

handles = list()
times = 0
for t in range(floor(server_threads / 10)):
    add = 0
    if db_all_rem > 0:
        add = 1
        db_all_rem -= 1

    my_range = db_all_mols[times * db_all_size:times * db_all_size + db_all_size + add]
   
    p = threading.Thread(target=measure, args=(my_range,))
    p.start()

    handles.append(p)
    times += 1

for handle in handles:
    handle.join()

