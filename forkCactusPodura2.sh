#!/bin/bash 
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=80G
#SBATCH --account=gompert-np
#SBATCH --qos=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=cactus
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

cd /scratch/general/nfs1/u6000989/cactus

module load cactus


max=2

jobs=(
#    "cactus jobStore5 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/cactusPodura_TpodE240140H1_TpodE240140H2.txt cactusPodura_TpodE240140H1_TpodE240140H2.hal --maxCores 20"
#    "cactus jobStore6 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/cactusPodura_TpodE240154H1_TpodE240154H2.txt cactusPodura_TpodE240154H1_TpodE240154H2.hal --maxCores 20"
#    "cactus jobStore7 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/cactusPodura_TpodE240140H1_TpodE240154H1.txt cactusPodura_TpodE240140H1_TpodE240154H1.hal --maxCores 20"
#    "cactus jobStore8 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/cactusPodura_TpodE240140H2_TpodE240154H2.txt cactusPodura_TpodE240140H2_TpodE240154H2.hal --maxCores 20"
#    "cactus jobStore9 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/cactusPodura_TpodE240140H1_TpodE240154H2.txt cactusPodura_TpodE240140H1_TpodE240154H2.hal --maxCores 20"
    "cactus jobStore10 /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/comp_aligns/cactusPodura_TpodE240140H2_TpodE240154H1.txt cactusPodura_TpodE240140H2_TpodE240154H1.hal --maxCores 20")

running=0

for job in "${jobs[@]}"; do
    eval "$job" &
    ((running++))

    # Once max jobs are running, wait for one to finish
    if (( running >= max )); then
        wait -n
        ((running--))
    fi
done

# Wait for any remaining jobs to finish
wait

