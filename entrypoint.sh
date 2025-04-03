#!/bin/bash

# This script directly executes the appropriate Python script based on the first argument
# No command name detection needed inside the container

# Debug information
echo "Running in container with arguments: $@"
echo "Working directory: $(pwd)"

# First argument should indicate which script to run
if [ "$1" = "micum" ]; then
    # Shift to remove the command name
    shift

    # Next argument is the input file path
    INPUT_FILE="$1"
    shift  # Remove the input file from arguments

    if [ -z "$INPUT_FILE" ]; then
        echo "Error: Input file not provided"
        echo "Usage: micum <input_file> [options]"
        exit 1
    fi

    if [ ! -f "$INPUT_FILE" ]; then
        echo "Error: Input file '$INPUT_FILE' not found"
        exit 1
    fi

    echo "Running MICUM.py with input file: $INPUT_FILE"
    cd /app
    python3 /app/MICUM.py "$INPUT_FILE" "$@"

elif [ "$1" = "merge_data" ]; then
    # Shift to remove the command name
    shift

    echo "Running merge_data.py with arguments: $@"
    cd /app
    python3 /app/merge_data.py "$@"

else
    echo "Unknown command: $1"
    echo "Supported commands: micum, merge_data"
    exit 1
fi