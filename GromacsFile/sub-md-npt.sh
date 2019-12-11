#!/bin/bash
#PBS -N n-ra
##PBS -M renh@upc.edu.cn
##PBS -m abe
##PBS -o /path/to/stdout
##PBS -e /path/to/stderr

#PBS -q bench
#PBS -l nodes=1:ppn=20
#PBS -l walltime=9:00:00

echo -n 'We work on:'
cat $PBS_NODEFILE
cd $PBS_O_WORKDIR

source /opt/gromacs/4.6.7/bin/GMXRC

#PYTHON=/opt/Anaconda/anaconda/bin/python
#SCRIPT=/home/ren/scripts/send_job_info.py


mdrun -nt 20 -pin on -deffnm npt 
#$PYTHON $SCRIPT -j $PBS_JOBNAME -I $PBS_JOBID -s terminates -d $PBS_O_WORKDIR
