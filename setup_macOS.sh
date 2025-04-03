#!/bin/bash

# MacOS Setup Script for MICUM

echo "Setting up MICUM for MacOS..."

# Create directory for scripts if it doesn't exist
SCRIPT_DIR="$HOME/bin"
mkdir -p "$SCRIPT_DIR"

# Add to PATH if not already there
if [[ ":$PATH:" != *":$SCRIPT_DIR:"* ]]; then
    echo "Adding $SCRIPT_DIR to PATH"
    echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.zshrc"
    echo 'export PATH="$HOME/bin:$PATH"' >> "$HOME/.bash_profile"
    export PATH="$HOME/bin:$PATH"
fi

# Create the micum script
cat > "$SCRIPT_DIR/micum" << 'EOF'
#!/bin/bash

# Check if we have arguments
if [ $# -eq 0 ]; then
    echo "Usage: $0 <input_file> [options]"
    exit 1
fi

# Get absolute path of the input file
INPUT_FILE=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
INPUT_DIR=$(dirname "$INPUT_FILE")
FILENAME=$(basename "$INPUT_FILE")

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found"
    exit 1
fi

# Remove first argument (input file) from args list
shift

# Run Docker container directly with the appropriate command
docker run --rm -it \
    -v "$INPUT_DIR:/workdir" \
    -w /workdir \
    micum \
    python3 /app/MICUM.py "/workdir/$FILENAME" "$@"
EOF

# Create the merge_data script
cat > "$SCRIPT_DIR/merge_data" << 'EOF'
#!/bin/bash

# Check if we have arguments
if [ $# -eq 0 ]; then
    echo "Error: No arguments provided"
    echo "Usage: $0 [options]"
    echo "Options should include -q and -m for input files"
    exit 1
fi

# Current directory for mounting
INPUT_DIR=$(pwd)

# Run Docker container directly with the appropriate command
docker run --rm -it \
    -v "$INPUT_DIR:/workdir" \
    -w /workdir \
    micum \
    python3 /app/merge_data.py "$@"
EOF

# Make scripts executable
chmod +x "$SCRIPT_DIR/micum"
chmod +x "$SCRIPT_DIR/merge_data"

echo "Building Docker image..."
docker build -t micum .

echo "Setup complete!"
echo "You can now use 'micum' and 'merge_data' commands."
echo "If commands are not found, please restart your terminal or run:"
echo "source ~/.zshrc  # for zsh"
echo "source ~/.bash_profile  # for bash"