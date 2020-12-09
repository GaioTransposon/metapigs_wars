
## usage example: 
## export input_file=/shared/homes/12705859/out_new_without_depth/contigs/contigs_withheaders/concat_assemblies_141.fa
## export out_dir=/shared/homes/12705859/macrel_all/out_contigs141
## qsub -V macrel_all.sh

#!/bin/bash
#PBS -l ncpus=20
#PBS -l walltime=120:00:00
#PBS -l mem=40g
#PBS -N macrel_all
#PBS -M daniela.gaio@student.uts.edu.au

which R 

source activate macrel_env

PATH="/shared/homes/12705859/miniconda3/envs/macrel_env/bin:$PATH"

which R

macrel contigs \
    --fasta $input_file \
    --output $out_dir \
    -t 18


# chmod 755 macrel_all.sh
# 
# export input_file=/shared/homes/12705859/out_new_without_depth/contigs/contigs_withheaders/concat_assemblies_141.fa
# export out_dir=/shared/homes/12705859/macrel_all/out_contigs141
# qsub -V macrel_all.sh
# 
# export input_file=/shared/homes/12705859/out_new_without_depth/contigs/contigs_withheaders/concat_assemblies_142.fa
# export out_dir=/shared/homes/12705859/macrel_all/out_contigs142
# qsub -V macrel_all.sh
# 
# export input_file=/shared/homes/12705859/out_new_without_depth/contigs/contigs_withheaders/concat_assemblies_143.fa
# export out_dir=/shared/homes/12705859/macrel_all/out_contigs143
# qsub -V macrel_all.sh
# 
# export input_file=/shared/homes/12705859/out_new_without_depth/contigs/contigs_withheaders/concat_assemblies_296.fa
# export out_dir=/shared/homes/12705859/macrel_all/out_contigs296
# qsub -V macrel_all.sh
# 
# export input_file=/shared/homes/12705859/out_new_without_depth/contigs/contigs_withheaders/concat_assemblies_297.fa
# export out_dir=/shared/homes/12705859/macrel_all/out_contigs297
# qsub -V macrel_all.sh
# 
# export input_file=/shared/homes/12705859/out_new_without_depth/contigs/contigs_withheaders/concat_assemblies_298.fa
# export out_dir=/shared/homes/12705859/macrel_all/out_contigs298
# qsub -V macrel_all.sh
# 
# export input_file=/shared/homes/12705859/out_new_without_depth/contigs/contigs_withheaders/concat_assemblies_299.fa
# export out_dir=/shared/homes/12705859/macrel_all/out_contigs299
# qsub -V macrel_all.sh