#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2024-2026 <actual copyright holder(s)>
# SPDX-License-Identifier: GPL-3.0-only
#
# This file is part of Phylo-MIP.
# See the LICENSE file in the project root for the full license text.

# Debug information
echo "Running in container"
echo "Working directory: $(pwd)"
echo "Arguments: $@"

# Ensure proper execution of Python scripts regardless of how they're called
if [[ "$1" == *"Phylo-MIP.py"* ]] || [[ "$2" == *"Phylo-MIP.py"* ]]; then
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