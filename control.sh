#!/bin/bash

# Define usage function
function usage {
  echo "Usage: $0 [command] [arguments]"
  echo "Options:"
  echo "  preprocess [arguments]    Run the preprocessing script with the provided arguments. Build a reference genome, count stuff, filter stuff with samtools."
  echo "  cellTag [arguments]    Run the cellTagR pipeline with the provided arguments"
  echo "  network [arguments]    Visualize clone lineages through networks and more."
  echo "  --help    Print this help. For more specific help on commands type [command] --help"
  exit 1
}

# Check if argument is provided
if [ $# -eq 0 ]; then
    echo "No Arguments provided"
    echo "Execute with --help for a list of arguments"
    exit 1
fi

# Parse arguments
case $1 in
  preprocess)
    echo "scripts/preprocessing.sh '${@:2}'" >> data/Log.log
    sh scripts/preprocessing.sh "${@:2}"
    ;;
  cellTag)
  echo "scripts/preprocessing.sh '${@:2}'" >> data/Log.log
    sh scripts/cellTag.sh "${@:2}"
    ;;
  network)
  echo "scripts/network.R '${@:2}'" >> data/Log.log
    sh scripts/network.R "${@:2}"
    ;;
  --help)
    usage
    ;;
  *)
    echo "Invalid Argument $1"
    echo "Execute with --help for a list of arguments"
    exit 1
    ;;
esac

