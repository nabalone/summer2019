#!/bin/bash
#
#SBATCH --job-name=gaus_weighttest
#SBATCH --output=gaussian_weighttest.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=1000:00
#SBATCH --mem-per-cpu=20000

cd /n/holylfs03/LABS/berger_lab/nchou 
python cnn_gaussianweightingtest.py


###END OF THE FILE#####


