# TpoduraColorGenetics
Research on the genetic basis of color in *T. podura* based on comparative genomics, GWA mapping and population genomic analyses

# Overview
My goal is to resolve the genetic basis of color pattern morphs in *T. podura*, with a focus on how this might differ at lower versus higher elevation sites and how it differs from (or not) what we see in *T. cristinae* and *T. chumahs*. Given what we have seen in this species and what we have already done in *T. podura*, my expectation is that chromosomal rearrangemetns (i.e., structural variation or SV) will be involved (associated or causal) (e.g., [Comeault et al. 2016](https://academic.oup.com/evolut/article-abstract/70/6/1283/6852167), [Villoutreix et al. 2020](https://www.science.org/doi/abs/10.1126/science.aaz4351)). My plans for this include GWA of color (using new and existing data), comparative alignments of multiple phased *T. podura* genomes and some local PCA analyses (among other things that remain TBD).

# Genomes and comparative alignments
We currenlty have three (x2) phased *T. podura* genomes, with more in the works. I don't think any of these three have been published yet. The first of these is from Dovetail (CEN4021), whereas all of the rest are (and will be) from Edinburgh. We expect to have at least 8x2 = 16 *T. podura genomes* (that is, at lest five more indiviudals that what we have at present). Everything is in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes` or `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edinburgh`. Here ist he breakdown of what we have:

| ID  | Location | Phenotype | Cactus | Annotation | 
|---------|-----|---------|:-:|:-:|
| cen4121 h1 | XX | Melanic | | |
| cen4121 h2 | XX | Melanic | | |
| 24_0140 h1 | DZR A | Melanic | | |
| 24_0140 h2 | DZR A | Melanic | | |
| 24_0154 h1 | BMT C | Green | | |
| 24_0154 h2 | BMT C | Green | | |

My first step with each genome is to split the fasta into files per haplotype and then to run repeat masking. This is done with `repeatmasker` (version 4.0.7); here is an example with 24_0140 and 24_0154:

```bash
#!/bin/sh 
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=384000
#SBATCH --account=gompert-np
#SBATCH --qos=gompert-np
#SBATCH --partition=gompert-np
#SBATCH --job-name=repeat
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zach.gompert@usu.edu

module load repeatmasker

#version 4.0.7
cd /uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/repeat_mask

## run repeat masker on each genome sequencei
## uses library from the 2020 Science paper developed by Victor

MAX_JOBS=4

genomes=(
	"/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0140/Hap1Chr.fasta"
	"/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0140/Hap2Chr.fasta"
	"/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0154/Hap1Chr.fasta"
	"/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edingburgh/24_0154/Hap2Chr.fasta"
)

for file in "${genomes[@]}"; do

	RepeatMasker -s -e ncbi -xsmall -pa 10 -lib RepeatLibMergeCentroidsRM.lib $file &
	# Limit the number of background jobs
	while (( $(jobs -rp | wc -l) >= MAX_JOBS )); do
		wait -n
  	done
done

# Wait for all remaining background jobs to finish
wait
```
