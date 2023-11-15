#! /bin/bash -login

## Exit if any command fails
set -e

## Load required modules
module load python/3.6.6

## Activate virtual environment
source env/bin/activate

case $1 in

    '-h' | '--help' | ' ' | '')
            echo -e '\e[31mSpecify which workflow to unlock (i.e. QTLtools_eQTL)'
            exit 2
            ;;
    
       'QTLtools_eQTL' | 'runQTLtools_eQTL')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/QTLtools_eQTL.smk --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status ./scripts/status.py
               ;;
        'responseQTL' | 'run_responseQTL')
                snakemake -j 1 --unlock -s snakefiles/responseQTL.smk --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status ./scripts/status.py   
            ;;

        'colocalization' | 'run_colocalization')
            snakemake -j 1 --unlock -s snakefiles/colocalization.smk --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status ./scripts/status.py
            ;;
esac

## Deactivate virtual environment
deactivate

## Success message
echo -e "\033[0;32mDirectory unlocked, ready to rerun."