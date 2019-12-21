#!/bin/bash
#
#SBATCH --job-name=b_50_mask__n_200_a_c_mp_5_kfold_0_l_0.0025
#SBATCH --output=cnn_run_b_50_mask__n_200_a_c_mp_5_kfold_0_l_0.0025.log
#SBATCH -p shared
#SBATCH --ntasks=1
#SBATCH --time=2000:00
#SBATCH --mem-per-cpu=10000

cd /n/home00/nchou/ia_work 
python cnn_aardvark_aug_concat.py -b 50 --mask  -n 200 -a -c --mp 5 --kfold 0 -l 0.0025 

###END OF THE FILE#####


