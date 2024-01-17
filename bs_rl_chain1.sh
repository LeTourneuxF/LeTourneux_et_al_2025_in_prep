#!/bin/bash
#SBATCH --time=167:00:00
#SBATCH --account=def-legagn3
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --output=chain1_p_het_Npheno_waic.out
#SBATCH --job-name=c1_Npheno_phet_waic
#SBATCH --mail-user=frederic.letourneux.1@ulaval.ca
#SBATCH --mail-type=ALL

module load gcc/9.3.0
module load r/4.2.1

Rscript joint_analysis_rldata_p_het_N_pheno_chain1.R