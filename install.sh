#!/bin/bash

# Download the micum script
curl -L https://raw.githubusercontent.com/yourusername/micum/main/micum -o /tmp/micum

# Make it executable
chmod +x /tmp/micum

# Move to a directory in PATH
sudo mv /tmp/micum /usr/local/bin/

# Pull the Docker image
docker pull shigebio/micum:latest

echo "MICUM installed successfully!"