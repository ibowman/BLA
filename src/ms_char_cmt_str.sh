#!/bin/bash

if [ ! -n "$3" ]
then
    echo "Usage: `basename $0` <wildcard> <injection-sites> <venv-dir>"
    exit 1
fi

files=$1
inj_sites=$2
venv_dir=$3

for file in $files
do
    name_base=`basename $file` # for output file
    logs_dir=qsub_logs
    job_name=cc-$name_base
    log_name_base=$logs_dir/${job_name}
    qsub -N $job_name -cwd -q compute.q -o $log_name_base.out -e $log_name_base.err -l hostslots=2,h_vmem=16G src/qsub_char_cmt_str.sh $file "$inj_sites" $venv_dir
done
