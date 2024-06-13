#!/bin/bash

queries=("DB00529.mol2" "DB00331.mol2" "DB01352.mol2" \
    "DB01621.mol2" "DB00728.mol2" "DB00887.mol2" \
"DB00050.mol2" "DB00403.mol2" "DB00732.mol2")

rm -rf ".jobs"; mkdir -p ".jobs"
rm -rf ".outputs"; mkdir -p ".outputs"

for query in ${queries[@]}; do
    job_file=".jobs/${query}.slurm"
    
    echo "#!/bin/bash
#SBATCH --job-name=pOptiPharm-${query}
#SBATCH --output=./outputs/${query}.txt

#SBATCH --partition=cpu_zen2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=96

# Script
hostname

cd /home/jero/project-optipharm/poptipharm/bin/comps
pwd

module purge
module load python py-pip gnu9 py-six

pip3 install numpy six pandas

python3.8 -u ./comps.py ${query}" > $job_file
    
    sbatch $job_file
done