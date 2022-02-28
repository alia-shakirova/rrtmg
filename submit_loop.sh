#!/bin/bash
project_name='062016_a10'
#output_scrt='/scratch/aliia/'${project_name}'/'
output_scrt='/home/aliia/tmp/new/'${project_name}'/'

if [ ! -d "${output_scrt}" ]; then
  mkdir -p "${output_scrt}"
fi

for t in {1..240}  
do
    if [ ! -d ${output_scrt}"$t" ]; then
      mkdir ${output_scrt}"$t"
    fi
    
    cd ${output_scrt}"$t"
    echo $(pwd)
    cp /home/aliia/RRTM/SW/script/albedo_feedback/main.m .
    cp /home/aliia/RRTM/SW/script/albedo_feedback/submit_job.pbs .
    sed -i "6c \ \ tind = $t;" main.m
    sed -i "5c #PBS -N ${project_name}_$t" submit_job.pbs
    sed -i "8c cd ${output_scrt}"$t"" submit_job.pbs
   # sed -i "7c  export PATH=$PATH:${output_scrt}"$t"" submit_job.pbs
    sed -i "11c matlab -nodisplay -nodesktop -singleCompThread < main.m > log_$t.out" submit_job.pbs
    qsub submit_job.pbs
done
