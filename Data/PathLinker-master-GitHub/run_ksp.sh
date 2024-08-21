#!/bin/bash

#SBATCH --job-name=PanCANPathLinkerBreastM2
#SBATCH --partition=gpu
#SBATCH --time=240:00:00
#SBATCH --gres=gpu:p100:4 
#SBATCH --ntasks=28
#SBATCH --ntasks-per-core=1
#SBATCH --mem=100g
#SBATCH --exclusive

# Check if the input file exists
if [ ! -f PanCandiffGenePathLinkerPairsUnique.txt ]; then
    echo "Error: File 'PanCandiffGenePathLinkerPairsUnique.txt' not found."
    exit 1
fi

# Read pairs from the file and store them in an array
mapfile -t pairs < PanCandiffGenePathLinkerPairsUnique.txt

# Loop through each pair
for pair in "${pairs[@]}"; do
    # Extract G and E from the pair
    G=$(echo "$pair" | awk '{print $1}')
    E=$(echo "$pair" | awk '{print $2}')

    # Define the output file name
    output_file="${G}_${E}_output.txt"

    # Run the command and redirect output to the output file
    python ksp_Astar.py PathLinkerGithubPPI.txt "$G" "$E" > "$output_file"
done

