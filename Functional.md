Run emapper.py (eggNOG-mapper v2)
1. Make sure your environment is activated
```bash
conda activate eggnogv2
```
create directory for eggNOG database properly
```bash
mkdir -p /home/bulat/miniconda3/lib/python3.12/site-packages/eggnogmapper/data
```
```bash
download_eggnog_data.py --data_dir /home/bulat/miniconda3/lib/python3.12/site-packages/eggnogmapper/data
```
uncompress the eggnog_proteins.dmnd.gz before emapper.py can be used.

Use of .faa file from gmhmmmp_*.faa file to run eggNOG-mapper for each sample! batch script
```bash
nano run_eggnog_mapper.sh
```
```bash
#!/bin/bash

# Base directory containing AA_* folders
base_dir="/mnt/e/Krona_results/below"

# Output directory for EggNOG results
output_dir="/mnt/e/Krona_results/eggnog_mapper_output"
mkdir -p "$output_dir"

# Path to EggNOG-mapper data directory (update if yours is different)
data_dir="/home/bulat/miniconda3/lib/python3.12/site-packages/eggnogmapper/data"


# Loop through all subfolders with FAA files
for folder in "$base_dir"/AA_*; do
    faa_file=$(find "$folder" -type f -name "*.faa")
    
    if [ -f "$faa_file" ]; then
        folder_name=$(basename "$folder")
        output_prefix="${output_dir}/${folder_name}_eggnog"

        echo "✅ Processing: $faa_file"
        
        emapper.py \
            -i "$faa_file" \
            --itype proteins \
            --output "$output_prefix" \
            --output_dir "$output_dir" \
            --data_dir "$data_dir" \
            --cpu 4
    else
        echo "❌ No FAA file found in $folder"
    fi
done
```
Using DIAMOND v2.1.6 to run emapper.py, run this script. 
output -- *_eggnog.emapper.annotations, *_eggnog.emapper.hits, *_eggnog.emapper.seed_orthologs
eggnog emapper annotations provides main functional annotation output. to be used for pathway-level analysis, GO terms, KEGG, eggNOG categories (COGs), etc
emapper.hits provides raw DIAMOND alignment (can be ignored unless troubleshooting)
emapper_seed_orthologs provide closest ortholog matches - to be used for fine-level comparative analysis

For the functional analysis, we would like to focus on:
1. KEGG orthology (KO)
2. eggNOG COG categories
3. Focus on GO terms
4. 


