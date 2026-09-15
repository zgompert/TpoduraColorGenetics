# TpoduraColorGenetics
Research on the genetic basis of color in *T. podura* based on comparative genomics, GWA mapping and population genomic analyses

# Overview
My goal is to resolve the genetic basis of color pattern morphs in *T. podura*, with a focus on how this might differ at lower versus higher elevation sites and how it differs from (or not) what we see in *T. cristinae* and *T. chumahs*. Given what we have seen in this species and what we have already done in *T. podura*, my expectation is that chromosomal rearrangemetns (i.e., structural variation or SV) will be involved (associated or causal) (e.g., [Comeault et al. 2016](https://academic.oup.com/evolut/article-abstract/70/6/1283/6852167), [Villoutreix et al. 2020](https://www.science.org/doi/abs/10.1126/science.aaz4351)). My plans for this include GWA of color (using new and existing data), comparative alignments of multiple phased *T. podura* genomes and some local PCA analyses (among other things that remain TBD). The main project is in `/uufs/chpc.utah.edu/common/home/gompert-group4/projects/tpodura_color_genomics`.

# Genomes and comparative alignments
We currenlty have three (x2) phased *T. podura* genomes, with more in the works. I don't think any of these three have been published yet. The first of these is from Dovetail (CEN4021), whereas all of the rest are (and will be) from Edinburgh. We expect to have at least 8x2 = 16 *T. podura genomes* (that is, at lest five more indiviudals that what we have at present). Everything is in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/t_podura_hap_cen4121` or `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/edinburgh`. Here ist he breakdown of what we have:

| ID  | Location | Phenotype | Cactus | Annotation | 
|---------|-----|---------|:-:|:-:|
| [cen4121 h1](https://github.com/user-attachments/files/32068438/CEN4121h1_report.html) | XX | Melanic | | |
| [cen4121 h2](https://github.com/user-attachments/files/32068436/CEN4121h2_report.html) | XX | Melanic | | |
| 24_0140 h1 | DZR A | Melanic | | |
| 24_0140 h2 | DZR A | Melanic | | |
| 24_0154 h1 | BMT C | Green | | |
| 24_0154 h2 | BMT C | Green | | |

My first step with each genome is to split the fasta into files per haplotype and then to run repeat masking. This is done with `repeatmasker` (version 4.0.7); here is an example with 24_0140 and 24_0154 (see `RunRMPodura26a.sh` in `/uufs/chpc.utah.edu/common/home/gompert-group4/data/timema/hic_genomes/repeat_mask`:

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
See also [forkCactusPodura.sh](forkCactusPodura.sh) and [forkCactusPodura2.sh](forkCactusPodura2.sh), as well as [batchHalPod.sh](batchHalPod.sh) and [HalForkPod.pl](HalForkPod.pl) for constructing the psl files (neede for dot plots) from the hal alignment files.

I moved the output from the comparative alignments (along with copies of scripts) to `/uufs/chpc.utah.edu/common/home/gompert-group4/projects/tpodura_color_genomics/comp_aligns`.

I then created synteny heatmaps and dot plots in R, see [SynPlotsPodura.R](SynPlotsPodura.R) (this includes some earlier comparisons with *T. cristinae*). My take thus far is that only h1 from 24_0154 is a ``green" allele, as it harbors a large deletion on chromosome 8 relative to the rest. But we will see if this holds with more data. Here are the dot plots I have so far: [AlnTpod_TpodH1_TpodE240140H1.pdf](https://github.com/user-attachments/files/32255231/AlnTpod_TpodH1_TpodE240140H1.pdf)
[AlnTpod_TpodH1_TpodE240154H1.pdf](https://github.com/user-attachments/files/32255230/AlnTpod_TpodH1_TpodE240154H1.pdf)
[AlnTpod_TpodH1_TpodE240140H2.pdf](https://github.com/user-attachments/files/32255228/AlnTpod_TpodH1_TpodE240140H2.pdf)
[AlnTpod_TpodH1_TpodE240154H2.pdf](https://github.com/user-attachments/files/32255227/AlnTpod_TpodH1_TpodE240154H2.pdf)
[AlnTpod_TpodE240140H1_podE240154H2.pdf](https://github.com/user-attachments/files/32255226/AlnTpod_TpodE240140H1_podE240154H2.pdf)
[AlnTpod_TpodE240140H2_podE240154H1.pdf](https://github.com/user-attachments/files/32255225/AlnTpod_TpodE240140H2_podE240154H1.pdf)
[AlnTpod_TpodE240140H2_podE240154H2.pdf](https://github.com/user-attachments/files/32255224/AlnTpod_TpodE240140H2_podE240154H2.pdf)
[AlnTpod_TpodE240154H1_podE240154H2.pdf](https://github.com/user-attachments/files/32255223/AlnTpod_TpodE240154H1_podE240154H2.pdf)
[AlnTpod_TpodE240140H1_podE240140H2.pdf](https://github.com/user-attachments/files/32255222/AlnTpod_TpodE240140H1_podE240140H2.pdf)
[AlnTpod_TpodE240140H1_podE240154H1.pdf](https://github.com/user-attachments/files/32255220/AlnTpod_TpodE240140H1_podE240154H1.pdf)


Next steps are:

1. GWA with existing GBS data, some combination of BSC, BMTC and IVC.
2. Sequence one moe green genome from low on the mountain.
3. Sequence four additional genomes from higher up on the mountain, including intermediate (yellow) and red morphs.
