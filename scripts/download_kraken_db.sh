#!/bin/bash

# Define the database links
16s_greengenes="https://genome-idx.s3.amazonaws.com/kraken/16S_Greengenes13.5_20200326.tgz"
16s_rdp="https://genome-idx.s3.amazonaws.com/kraken/16S_RDP11.5_20200326.tgz"
16s_silva="https://genome-idx.s3.amazonaws.com/kraken/16S_Silva132_20200326.tgz"


# Parse the arguments
for i in "$@"
do
case $i in
    --method=*)
    method="${i#*=}"
    shift
    ;;
    --references=*)
    references="${i#*=}"
    shift
    ;;
    --output_dir=*)
    output_dir="${i#*=}"
    shift
    ;;
    *)
    ;;
esac
done

# Select the appropriate method

if [ "$method" == "target" ]; then

  if [ "$references" == "greengenes" ]; then
      link=${16s_greengenes}
      wget -O ${link} | gunzip | tar -xvf -C ${output_dir}
  elif [ "$references" == "silva" ]; then
      link=${16s_silva}
      wget -O ${link} | gunzip | tar -xvf -C ${output_dir}
  elif [ "$references" == "rdp" ]; then
      link=${16s_rdp}
      wget -O ${link} | gunzip | tar -xvf -C ${output_dir}
  else
    echo "Invalid reference argument. Please specify either 'greengenes' or 'rdp' or silva"
    exit 1
  fi
elif [ "$method" == "shotgun" ]; then
  kraken2-build --download-taxnomy --db ${output_dir}
  kraken2-build --download-library ${references} --db ${output_dir}
  kraken2-build --build --db ${output_dir}
else
  echo "invalid method argument. Please specify either 'target' or 'shotgun'"
  exit 1
fi





