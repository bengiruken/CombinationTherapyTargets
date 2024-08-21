#!/bin/bash

#SBATCH --job-name=EnrichRPanCan
#SBATCH --partition=gpu
#SBATCH --time=240:00:00
#SBATCH --gres=gpu:p100:4 
#SBATCH --ntasks=28
#SBATCH --ntasks-per-core=1
#SBATCH --mem=100g
#SBATCH --exclusive


python EnrichR_Shortest_Path_Biowulf.py


