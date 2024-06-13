#!/bin/bash

#SBATCH -J pOptiData
#SBACTH -D /home/jero/project-optipharm/poptipharm/bin

#SBATCH -o /home/jero/project-optipharm/poptipharm/bin/data.output

#SBATCH -p cpu_zen2
#SBATCH -N 1
#SBATCH --cpus-per-task=96
#SBATCH --nodelist=c25

hostname

cd /home/jero/project-optipharm/poptipharm/bin
pwd

module purge
module load python py-pip gnu9 py-six

pip3 install numpy six pandas

python3.8 -u ./data.py