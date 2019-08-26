#!/bin/bash
#
#SBATCH --job-name=cnn_sample_compy
#SBATCH --output=cnn_sample_compy.log
#
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=1000:00
#SBATCH --mem-per-cpu=10000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_sample_compy.py


###END OF THE FILE#####


