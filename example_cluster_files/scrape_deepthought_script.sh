#!/bin/bash

#SBATCH -t 120:00
#SBATCH --mem=15000
#SBATCH --ntasks=1
#SBATCH --share


cd /lustre/asavel/wasp_76

source wasp_76_env/bin/activate
python scrape_deepthought_data.py
