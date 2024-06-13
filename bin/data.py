import subprocess
import shutil
import pandas as pd
import statistics

db_mol = ["DB00460", "DB00093", "DB00014"]
db_threads = [1, 2, 4, 8, 16, 32, 64, 95]
db_specs = [4, 8, 16, 32, 64, 128]
db_levels = [5, 10, 15, 20, 25]
db_ratio = [1.0, 2.0, 4.0]

tries = 10

initial_evals = 200000

specs_evals_ratio = initial_evals / db_specs[0]

header = ["OP Version", "Query", "Target", "Threads", "Pinned", "Max. Evals", "Max. Specs.", "Levels", "Initial angles", "Avg. Value", "Std. Dev. Value", "Avg. Time (s)", "Std. Dev. Time", "Speed-up relative secuential"]
table = pd.DataFrame(columns=header)

shutil.copy("./OPini.FMOL", "./OPini.JERO")

edit_cmd = "./pOPShapeSimilarity -c JERO -N {} -M {} -L {}"

cmd = "./OPShapeSimilarity -c JERO -q ./../../DB/fda/{}.mol2 -d ./../../DB/fda/{}.mol2 -h 1 -S"
parallel_cmd = "./pOPShapeSimilarity -c JERO -q ./../../DB/fda/{}.mol2 -d ./../../DB/fda/{}.mol2 -h 1 -S -th {} -pn {}"

for query in db_mol:
    for target in db_mol:
        if query == target:
            continue

        print("Comparing {} vs {}.".format(query, target))

        for th in db_threads:
            print("\tCalculating for {} threads.".format(th))
            for pinned in range(2):
                print("\t\tPinned {}.".format(pinned))
                for evals_ratio in db_ratio[:th+1]:
                    print("\t\tNew evals. ratio: {}.".format(evals_ratio))
                    for spec_num in db_specs[:th+1]:
                        print("\t\tMax. specs: {}.".format(spec_num))
                        for lvl in db_levels[:th+1]:
                            print("\t\tMax. level: {}.".format(lvl))
                            current_evals = specs_evals_ratio * spec_num * evals_ratio
                            opini = subprocess.Popen(edit_cmd.format(current_evals, spec_num, lvl), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, universal_newlines=True) # Changed OPini.JERO
                            opini, err_opini = opini.communicate()

                            fmt_sec_cmd = cmd.format(query, target)
                            fmt_par_cmd = parallel_cmd.format(query, target, th, pinned)

                            sec_values = list()
                            sec_times = list()

                            par_values = list()
                            par_times = list()

                            for i in range(tries):
                                print("\t\t\tTry {}.".format(i + 1))
                                if pinned == 0 and th < 2:
                                    sec_output = subprocess.Popen(fmt_sec_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, universal_newlines=True)
                                    
                                    sec_output, err_a = sec_output.communicate()
                                    
                                    # Parse secuential output
                                    sec_lines = sec_output.splitlines()
                                    # Sec data on before last row
                                    sec_data = sec_lines[0].split()

                                    # Get time (pos 5)
                                    sec_time = sec_data[4]

                                    # Get value (pos 17)
                                    sec_value = sec_data[16]

                                    # Append
                                    sec_values.append(float(sec_value))
                                    sec_times.append(float(sec_time))

                                par_output = subprocess.Popen(fmt_par_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, universal_newlines=True)

                                par_output, err_b = par_output.communicate()

                                # Parse parallel output
                                par_lines = par_output.splitlines()
                                # Par data on before-before last row
                                par_data = par_lines[len(par_lines) - 2]
                                par_data = par_data.split()

                                # Get time (pos 5)
                                par_time = par_data[4]
                                
                                # Get value (pos 17)
                                par_value = par_data[16]

                                # Append
                                par_values.append(float(par_value))
                                par_times.append(float(par_time))

                            # New row
                            if pinned == 0 and th < 2:
                                table.loc[len(table.index)] = ["Secuential", query, target, 1, "-", current_evals, spec_num, lvl, max(1, min(12, spec_num)), statistics.mean(sec_values), statistics.stdev(sec_values), statistics.mean(sec_times), statistics.stdev(sec_times), "-"]
                                table.loc[len(table.index)] = ["Threads", query, target, th, pinned, current_evals, spec_num, lvl, max(1, min(12, spec_num)), statistics.mean(par_values), statistics.stdev(par_values), statistics.mean(par_times), statistics.stdev(par_times), statistics.mean(sec_times) / statistics.mean(par_times)]
                            else:
                                table.loc[len(table.index)] = ["Threads", query, target, th, pinned, current_evals, spec_num, lvl, max(1, min(12, spec_num)), statistics.mean(par_values), statistics.stdev(par_values), statistics.mean(par_times), statistics.stdev(par_times), "-"]
                            table.to_csv("./temp-data.csv", index=False)

table.to_csv("./all-data.csv", index=False)