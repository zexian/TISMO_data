#!/bin/bash
#SBATCH --job-name=geo1
#SBATCH --mem=10G       # total memory need
#SBATCH --time=96:00:00
#SBATCH -c 2 #number of cores
#SBATCH -o o_%j.txt
#SBATCH -e e_%j.txt


python scrna_parser_runner.py parser -d 2019/01/23-2019/05/23 -el ./problem_GSE.txt -o syng_collection2.xls -p geo_xml2  