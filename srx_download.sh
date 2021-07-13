#!/bin/bash
#SBATCH --job-name=vv_2dl
#SBATCH --mem=30G       # total memory need
#SBATCH --time=24:00:00
#SBATCH -c 2 #number of cores
#SBATCH -o dl_vv_o_%j.txt
#SBATCH -e dl_vv_e_%j.txt

python srx_download.py
