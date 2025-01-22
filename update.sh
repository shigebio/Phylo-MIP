#!/bin/bash

# If an error occurs, exit
set -e

# Project name and fixed image name
PROJECT_NAME="micum_project"
IMAGE_NAME="micum:latest"

# Stop and remove the old container
echo "Stopping and removing the existing containers..."
sudo docker-compose -p "$PROJECT_NAME" down

# Delete the old image
if sudo docker images "$IMAGE_NAME" | grep -q "$IMAGE_NAME"; then
  echo "Removing the old Docker image: $IMAGE_NAME"
  sudo docker rmi "$IMAGE_NAME"
fi

# Update the GitHub repository
echo "Updating the project from GitHub..."
git fetch origin
git reset --hard origin/main

# Build a new image using Docker Compose
echo "Building the Docker image..."
sudo docker-compose -p "$PROJECT_NAME" build

# Start a new container
echo "Starting the new containers..."
sudo docker-compose -p "$PROJECT_NAME" up -d

# Clean up unwanted images
echo "Cleaning up dangling images..."
sudo docker image prune -f

echo "Update completed for project: $PROJECT_NAME!"
