#!/bin/sh 
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=240G
#SBATCH --account=gompert
#SBATCH --qos=gompert-grn
#SBATCH --partition=gompert-grn
#SBATCH --job-name=cactus-syn
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu


cd /scratch/general/nfs1/u6000989/cactus

perl /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/HalForkPod.pl cactusPodura_TpodE24*hal
