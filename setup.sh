#!/bin/bash

# Simple setup script using Makefile
# Usage: ./setup.sh

set -e

echo "🚀 Setting up Phylo-MIP..."
echo ""

# Check if make is available
if ! command -v make &> /dev/null; then
    echo "❌ Error: 'make' is not installed"
    echo ""
    case "$(uname -s)" in
        Darwin)
            echo "On macOS, install with:"
            echo "  xcode-select --install"
            echo "  or install via Homebrew: brew install make"
            ;;
        Linux)
            echo "On Linux, install with:"
            echo "  sudo apt-get install build-essential  # Ubuntu/Debian"
            echo "  sudo yum install make                 # CentOS/RHEL"
            ;;
        *)
            echo "Please install 'make' for your system"
            ;;
    esac
    exit 1
fi

# Check if docker is available
if ! command -v docker &> /dev/null; then
    echo "❌ Error: Docker is not installed"
    echo ""
    echo "Please install Docker from: https://docs.docker.com/get-docker/"
    exit 1
fi

# Check if Docker is running
if ! docker info &> /dev/null; then
    echo "❌ Error: Docker is not running"
    echo "Please start Docker and try again."
    exit 1
fi

# Run make setup
make setup

echo ""
echo "🎉 Setup complete!"
echo ""
echo "Quick start:"
echo "  make help                              # Show all available commands"
echo "  make phylo-mip INPUT=your_file.fasta  # Run Phylo-MIP analysis"
echo ""