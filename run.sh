#!/bin/bash -e

TAXIDS=$1
FOLDER=$2
COVERAGE=$3
WINDOW_SIZE=$4
KMER=$5
THREADS=$6

GENMAP="./build/bin/genmap"
CLASSIFY="./build/bin/classify"
MCPD="../pipeline/decoder/MCPD"

READ_FOLDER="$FOLDER/reads_${COVERAGE}x"
DESIGN_FILE="${FOLDER}/genmap_${WINDOW_SIZE}_${KMER}.design"
KMER_FILE="${FOLDER}/genmap_${WINDOW_SIZE}_${KMER}.kmers"

#if [ -f "$FOLDER/taxids" ]; then
#  PREV_TAXIDS=$(cat "$FOLDER/taxids")
#  if [ "$PREV_TAXIDS" != "$TAXIDS" ]; then
#    echo "Running with different taxids than before in the same directory. Deleting everything and starting from scratch."
#    rm -rf "$FOLDER"
#  fi
#fi

function check_command {
  if ! [ -x "$(command -v $1)" ]; then
    echo "Error: $1 is not installed." >&2
    exit 1
  fi
}

check_command "bc"
check_command "wgsim"
check_command "parallel"

# Downloading the genomes
IFS=',' read -r -a taxids_array <<< "$TAXIDS"
for id in "${taxids_array[@]}"; do
  krakenuniq-download --db ${FOLDER} refseq/bacteria/Complete_Genome/species_taxid=${id} --threads ${THREADS}
done

# Indexing the genomes
if [ ! -d "$FOLDER/index" ]; then
  ${GENMAP} index -FD ${FOLDER} -I "${FOLDER}/index"
fi

# Computing the design file
if [ ! -f "${DESIGN_FILE}" ]; then
  ${GENMAP} map -I "${FOLDER}/index" --design --exclude-pseudo -fs -O "${FOLDER}" --design-window ${WINDOW_SIZE} -E 0 -K ${KMER} -T ${THREADS}
  mv "${FOLDER}/genmap.design" "${DESIGN_FILE}"
  mv "${FOLDER}/genmap.kmers" "${KMER_FILE}"
fi

# Simulating reads
if [ ! -d "$READ_FOLDER" ]; then
  echo "Simulating 150bp reads"

  mkdir -p "$READ_FOLDER"

  for d in $(find "${FOLDER}" -iname "*.fna"); do
    BASES=$(awk '!($0 ~ /^>/) { printf "%s", $0 }' $d |  wc -m)
    READS_F=$(echo "($COVERAGE + $BASES) / 150.0" | bc)
    READS=${READS_F%.*}
    filename=$(basename -- "$d")
    filename_wo_ext="${filename%.*}"
    
    #echo "$BASES $READS $filename"
    wgsim -N ${READS} -1 150 -2 150 $d "${READ_FOLDER}/${filename_wo_ext}.fq" /tmp/2nd_paired.fq > /dev/null 2>&1
    
    # fastq to fasta
    awk '(NR%4==1) { print ">" substr($0, 2) } (NR%4==2)' "${READ_FOLDER}/${filename_wo_ext}.fq" > "${READ_FOLDER}/${filename_wo_ext}.fa"
    rm /tmp/2nd_paired.fq "${READ_FOLDER}/${filename_wo_ext}.fq"
  done
fi

# Classify and decode reads
#for d in $(find "${READ_FOLDER}" -iname '*.fa'); do
#  ${CLASSIFY} -K ${KMER_FILE} -F "${d}"
#  ${MCPD} -w 1000 -s 10000 -m 0.001000 ${DESIGN_FILE} "${d}.result" > "${d}.decoded"
#done
find "${READ_FOLDER}" -iname '*.fa' | parallel --will-cite -j ${THREADS} "${CLASSIFY} -K ${KMER_FILE} -F {} && ${MCPD} -w 1000 -s 10000 -m 0.001000 ${DESIGN_FILE} {}.result > {}.decoded"

# Summarize results
find "${READ_FOLDER}" -iname '*.decoded' | xargs awk '{ if (last_filename == FILENAME) { if (!($0 ~ /^#/) && $2 > 0.02) { print } } else { if (NR!=1) { print "" } print FILENAME; last_filename = FILENAME } }' > "${FOLDER}/summary_${WINDOW_SIZE}_${KMER}"

# Print results
cat "${FOLDER}/summary_${WINDOW_SIZE}_${KMER}" | column -t -s $'\t'
