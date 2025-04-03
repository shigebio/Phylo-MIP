#!/bin/bash

# Simple entrypoint that selects the appropriate Python script based on the first argument

# Debug information
echo "Running in container"
echo "Working directory: $(pwd)"
echo "Arguments: $@"

# First argument determines the command to run
if [[ "$1" == *"MICUM.py"* ]]; then
    # Running the MICUM.py script
    echo "Executing MICUM.py with arguments: $@"
    python3 "$@"
elif [[ "$1" == *"merge_data.py"* ]]; then
    # Running the merge_data.py script
    echo "Executing merge_data.py with arguments: $@"
    python3 "$@"
else
    # Determine if we're being run with a direct command name
    if [ "$1" = "micum" ]; then
        shift
        INPUT_FILE="$1"
        shift
        echo "Executing MICUM.py with input file: $INPUT_FILE and options: $@"
        python3 /app/MICUM.py "$INPUT_FILE" "$@"
    elif [ "$1" = "merge_data" ]; then
        shift
        echo "Executing merge_data.py with options: $@"
        python3 /app/merge_data.py "$@"
    else
        echo "Unknown command or script: $1"
        echo "Available commands: micum, merge_data"
        echo "Or specify full path to script: /app/MICUM.py, /app/merge_data.py"
        exit 1
    fi
fi