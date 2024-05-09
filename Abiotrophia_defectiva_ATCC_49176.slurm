#!/bin/bash
#SBATCH --job-name=hostname
#SBATCH --partition=qcpu_23i   
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G  #
#SBATCH --mail-type=end
#SBATCH --mail-user=dengxg@tib.cas.cn
#SBATCH --output=/hpcfs/fhome/dengxg/virtual_metabolic_human/%j.out
#SBATCH --error=/hpcfs/fhome/dengxg/virtual_metabolic_human/%j.err

conda activate dxbase

cd /hpcfs/fproject/dengxg/mqc
python start3.py -m '/hpcfs/fproject/dengxg/mqc/mqc/local_test_data/virtual_metabolic_human/Abiotrophia_defectiva_ATCC_49176.xml' -o 'hpcfs/fproject/dengxg/mqc/tmp/virtual_metabolic_human/Abiotrophia_defectiva_ATCC_49176'