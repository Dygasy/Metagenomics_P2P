# Metagenomics_V2

Meganize a .daa file on ubuntu, you'll need MEGAN command-line tools, specifically daa-meganizer, which annotates your .daa file for MEGAN. From here, you will be able to: 
1. Access GTDB classification down to species (if annotated properly)
2. See functional classifications: EC, eggNOG, SEED, InterPro2GO
3. Export and analyse downstream in MEGAN and R
4. 
Install MEGAN 6.25.10: 
```bash
wget https://software-ab.informatik.uni-tuebingen.de/download/megan6/MEGAN_Community_unix_6_24_22.sh
bash MEGAN_Community_unix_6_24_22.sh
```
Download NCBI-based mapping file megan-map-Feb2022.db
Files are in E:/NHCS/NHCSMS/daa/megan
Create the script
In Ubuntu, run:

bash
Copy
Edit
nano batch_meganize.sh
Paste this:

```bash
#!/bin/bash

# ===== CONFIGURATION =====
DAA2RMA6="/mnt/c/nr_db/megan/tools/daa2rma"
MAP_DB="/mnt/c/nr_db/megan/megan-map-Feb2022.db/megan-map-Feb2022.db"
THREADS=36

BATCH_FOLDERS=(
    "/home/bulat/Dayang/batch_1"
    "/home/bulat/Dayang/batch_2"
)

convert_batch() {
    local folder="$1"
    echo "🟡 Processing: $folder"
    for daa in "$folder"/*.daa; do
        out="${daa%.daa}.rma6"
        echo "🔄 Converting: $daa"
        "$DAA2RMA6" \
            --in "$daa" \
            --out "$out" \
            --acc2taxa "$MAP_DB" \
            --mapDB "$MAP_DB" \
            --threads "$THREADS"
        echo "✅ Finished: $out"
    done
}

for folder in "${BATCH_FOLDERS[@]}"; do
    convert_batch "$folder" &
done

wait
echo "🏁 All RMA6 conversions done!"

```
2. Make it executable
```bash
chmod +x batch_meganize.sh
```
3. Run the script
```bash
./batch_meganize.sh
```
Workstation is able to process using 36-core, 72-thread Intel Xeon w9-3475X with 117GB RAM. 
