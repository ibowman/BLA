#!/bin/bash

if [ ! -n "$7" ]
then
    echo 'Usage: `basename $0` <char-cmt-str-csv> <case-dir> <grid-cell-size> <channel> "<inj-site-colors>" "<levels>" <venv-dir>'
    exit 1
fi

char_cmt_str_csv=$1;
case_dir=$2;
grid_cell_size=$3;
channel=$4;
inj_site_colors=$5
levels=$6
venv_dir=$7;

logs_dir=qsub_logs
mkdir -p $logs_dir # prepare logs file if it doesn't already exist
name_base=`basename $case_dir` # for output file

for lvl in $levels; do
    print_lvl=`printf "%0.3d" $lvl`
    job_name=cc${print_lvl}-${channel}-${name_base}
    log_name_base=$logs_dir/$job_name
    qsub -l hostslots=2,h_vmem=16G -N $job_name -cwd -o $log_name_base.out -e $log_name_base.err -q compute.q src/qsub_run_cmt_clr_thresh.sh $char_cmt_str_csv $case_dir $grid_cell_size $channel $lvl "$inj_site_colors" $venv_dir
done;
