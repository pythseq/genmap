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

ERROR_RATE="0.02" # for simulated reads

READ_FOLDER="$FOLDER/reads_${COVERAGE}x"
DESIGN_FILE="${FOLDER}/genmap_${WINDOW_SIZE}_${KMER}.design"
DESIGN_FILE_MODIFIED="${FOLDER}/genmap_${WINDOW_SIZE}_${KMER}.design.tmp"
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
    READS_F=$(echo "($COVERAGE * $BASES) / 150.0" | bc)
    READS=${READS_F%.*}
    filename=$(basename -- "$d")
    filename_wo_ext="${filename%.*}"
    
    #echo "$BASES $READS $filename"
    wgsim -N ${READS} -1 150 -2 150 -e $ERROR_RATE $d "${READ_FOLDER}/${filename_wo_ext}.fq" /tmp/2nd_paired.fq > /dev/null 2>&1
    
    # fastq to fasta
    awk '(NR%4==1) { print ">" substr($0, 2) } (NR%4==2)' "${READ_FOLDER}/${filename_wo_ext}.fq" > "${READ_FOLDER}/${filename_wo_ext}.fa"
    rm /tmp/2nd_paired.fq "${READ_FOLDER}/${filename_wo_ext}.fq"
  done
fi

# Compute FN rate
FN_RATE=$(echo "$ERROR_RATE $KMER $COVERAGE" | awk '{ printf "%.2f", (1 - ((1 - $1)^($2)))^($3) }')
TN_RATE=$(echo "$FN_RATE" | awk '{ printf "%.2f", (1 - $1) }')
echo "FN rate: ${FN_RATE}"

head -n 3 ${DESIGN_FILE} > ${DESIGN_FILE_MODIFIED}
echo -e "1.0\t${FN_RATE}\t${TN_RATE}" >> ${DESIGN_FILE_MODIFIED}
tail -n +5 ${DESIGN_FILE} >> ${DESIGN_FILE_MODIFIED}

# Classify and decode reads
# TODO: replace DESIGN_FILE with DESIGN_FILE_MODIFIED for computed FN rate (in the line below)
find "${READ_FOLDER}" -iname '*.fa' | parallel --will-cite -j ${THREADS} "${CLASSIFY} -K ${KMER_FILE} -F {} && ${MCPD} -w 1000 -s 10000 -m 0.001 ${DESIGN_FILE_MODIFIED} {}.result > {}.decoded"

# Summarize decoding results
find "${READ_FOLDER}" -iname '*.decoded' | xargs awk -F$'\t' -f "./decodingSummary.awk" > "${FOLDER}/summary_${WINDOW_SIZE}_${KMER}"

# Print coding results
cat "${FOLDER}/summary_${WINDOW_SIZE}_${KMER}" | column -t -s $'\t'

# Build kraken database
if [ ! -f "${FOLDER}/database.kdb" ]; then
  echo "Build kraken database"
  krakenuniq-build --build --db ${FOLDER} --threads ${THREADS}
fi

# Run krakenuniq
echo "Running krakenuniq (if not already run)"
for d in $(find "${READ_FOLDER}" -iname "*.fa"); do
  if [ ! -f "${d}.krakenuniq.report" ]; then
    krakenuniq --report-file ${d}.krakenuniq.report ${d} --db ${FOLDER} --threads ${THREADS} --fasta-input > /dev/null 2>&1
  fi
done


