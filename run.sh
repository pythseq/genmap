#!/bin/bash -xe

TAXIDS=$1
FOLDER=$2
THREADS=$3

GENMAP="./build/bin/genmap"

#if [ -f "$FOLDER/taxids" ]; then
#  PREV_TAXIDS=$(cat "$FOLDER/taxids")
#  if [ "$PREV_TAXIDS" != "$TAXIDS" ]; then
#    echo "Running with different taxids than before in the same directory. Deleting everything and starting from scratch."
#    rm -rf "$FOLDER"
#  fi
#fi

#mkdir -p "$FOLDER"

IFS=',' read -r -a taxids_array <<< "$TAXIDS"

for id in "${taxids_array[@]}"; do
  krakenuniq-download --db ${FOLDER} refseq/bacteria/Complete_Genome/species_taxid=${id} --threads ${THREADS}
done

if [ ! -d "$FOLDER/index" ]; then
  ${GENMAP} index -FD ${FOLDER} -I "${FOLDER}/index"
fi


