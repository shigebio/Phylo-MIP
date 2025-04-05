#!/bin/bash

# Debug information
echo "Running in container"
echo "Working directory: $(pwd)"
echo "Arguments: $@"

# Ensure proper execution of Python scripts regardless of how they're called
if [[ "$1" == *"MICUM.py"* ]] || [[ "$2" == *"MICUM.py"* ]]; then
    # Make sure we're using python3 to execute the script
    if [[ "$1" == "python3" ]]; then
        echo "Executing: $@"
        exec "$@"
    else
        echo "Executing: python3 $@"
        exec python3 "$@"
    fi
elif [[ "$1" == *"merge_data.py"* ]] || [[ "$2" == *"merge_data.py"* ]]; then
    # Make sure we're using python3 to execute the script
    if [[ "$1" == "python3" ]]; then
        echo "Executing: $@"
        exec "$@"
    else
        echo "Executing: python3 $@"
        exec python3 "$@"
    fi
else
    # Others
    echo "Executing command: $@"
    exec "$@"
fi