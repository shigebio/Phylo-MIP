#!/bin/bash

# Detect shell configuration file
if [[ "$SHELL" == */zsh ]]; then
  CONFIG_FILE="$HOME/.zshrc"
else
  CONFIG_FILE="$HOME/.bashrc"
fi

# Directory of the MICUM project (where docker-compose.yml is located)
MICUM_DIR="$(pwd)"  # Use the current directory as the MICUM project directory

# Check if Docker can run without sudo (i.e., user is in the docker group)
if docker ps >/dev/null 2>&1; then
  if command -v docker-compose >/dev/null 2>&1; then
    DOCKER_COM_CMD="docker-compose"
  else
    DOCKER_COM_CMD="docker compose"
  fi
else
  if command -v docker-compose >/dev/null 2>&1; then
    DOCKER_COM_CMD="sudo docker-compose"
  else
    DOCKER_COM_CMD="sudo docker compose"
  fi
fi

# Check if the Docker image already exists, and only build if necessary
if ! docker images | grep -q "micum"; then
  echo "Building Docker image..."
  (cd "$MICUM_DIR" && $DOCKER_COM_CMD build)
else
  echo "Docker image already exists. Skipping build."
fi

# Remove existing micum() function definition from shell config if it exists
if grep -q "micum()" "$CONFIG_FILE"; then
  echo "Removing existing micum() function..."
  if [[ "$(uname)" == "Darwin" ]]; then
    sed -i '' '/micum()/,/^}/d' "$CONFIG_FILE"
  else
    sed -i '/micum()/,/^}/d' "$CONFIG_FILE"
  fi
fi

# Add micum function to shell config (only if it hasn't been added yet)
if ! grep -q "micum()" "$CONFIG_FILE"; then
  echo "Adding micum function..."
  cat << 'EOF' >> "$CONFIG_FILE"

# MICUM command function
micum() {
    set -e
    set -x

    if docker ps >/dev/null 2>&1; then
      DOCKER_CMD="docker"
    else
      DOCKER_CMD="sudo docker"
    fi

    # Set the MICUM_DIR as an environment variable
    export MICUM_DIR="$(pwd)"

    if ! $DOCKER_CMD ps | grep -q micum; then
      echo "Starting MICUM container..."
      $DOCKER_CMD compose up -d
    fi

    echo "Running MICUM.py with arguments: $@"
    if ! $DOCKER_CMD exec -i -w /app micum env PYTHONUNBUFFERED=1 MICUM_DIR="$MICUM_DIR" python3 MICUM.py "$@"; then
      echo "Error: MICUM.py execution failed."
    fi
}
EOF
fi

# Apply changes
exec $SHELL
