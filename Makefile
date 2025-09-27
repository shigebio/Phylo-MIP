# Phylo-MIP Makefile
# Cross-platform Docker-based pipeline

# Detect platform for Docker
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    PLATFORM_FLAG = --platform linux/amd64
else
    PLATFORM_FLAG =
endif

# Detect user mapping (for file permissions)
USER_FLAG := $(shell if command -v id >/dev/null 2>&1; then echo "--user $$(id -u):$$(id -g)"; fi)

# Docker image name
IMAGE_NAME = phylo-mip

.PHONY: help setup build clean phylo-mip merge-data update

help:
	@echo "Phylo-MIP Pipeline Commands:"
	@echo ""
	@echo "Setup and maintenance:"
	@echo "  make setup     - Build Docker image and setup environment"
	@echo "  make build     - Build/rebuild Docker image"
	@echo "  make update    - Update and rebuild (alias for build)"
	@echo "  make clean     - Remove Docker image"
	@echo ""
	@echo "Usage:"
	@echo "  make phylo-mip INPUT=<file> [ARGS='additional args']"
	@echo "  make merge-data QIIME=<file> PHYLOMIP=<file> FORMAT=<format> [OUTPUT=<file>] [DEBUG=1]"
	@echo ""
	@echo "Examples:"
	@echo "  make phylo-mip INPUT=data.fasta"
	@echo "  make phylo-mip INPUT=data.fasta ARGS='--threads 4 --method fast'"
	@echo "  make merge-data QIIME=qiime_output.tsv PHYLOMIP=phylo_results.csv FORMAT=csv"
	@echo "  make merge-data QIIME=qiime.tsv PHYLOMIP=phylo.csv FORMAT=tsv OUTPUT=merged.tsv DEBUG=1"
	@echo ""
	@echo "Notes:"
	@echo "  - All file paths should be relative to current directory"
	@echo "  - Output files will be created in the same directory as input files"
	@echo "  - Use DEBUG=1 to enable debug mode for merge-data"

setup: build
	@echo ""
	@echo "✅ Setup complete!"
	@echo ""
	@echo "You can now use:"
	@echo "  make phylo-mip INPUT=<your_file>"
	@echo "  make merge-data QIIME=<qiime_file> PHYLOMIP=<phylomip_file> FORMAT=<format>"
	@echo ""
	@echo "Run 'make help' anytime to see all available commands."

build:
	@echo "🔨 Building Docker image..."
	@docker build $(PLATFORM_FLAG) -t $(IMAGE_NAME) .
	@echo "✅ Docker image built successfully!"

update: build
	@echo "✅ Update complete!"

clean:
	@echo "🧹 Removing Docker image..."
	@docker rmi $(IMAGE_NAME) 2>/dev/null || echo "Image already removed or doesn't exist"
	@echo "✅ Cleanup complete!"

# Phylo-MIP command
phylo-mip:
	@if [ -z "$(INPUT)" ]; then \
		echo "❌ Error: INPUT parameter is required"; \
		echo "Usage: make phylo-mip INPUT=<file> [ARGS='additional args']"; \
		echo "Example: make phylo-mip INPUT=data.fasta"; \
		exit 1; \
	fi
	@if [ ! -f "$(INPUT)" ]; then \
		echo "❌ Error: Input file '$(INPUT)' not found"; \
		echo "Please check the file path and try again."; \
		exit 1; \
	fi
	@echo "🧬 Running Phylo-MIP on $(INPUT)..."
	@echo "📁 Working directory: $$(cd "$$(dirname "$(INPUT)")" && pwd)"
	@docker run --rm -i $(PLATFORM_FLAG) $(USER_FLAG) \
		-v "$$(cd "$$(dirname "$(INPUT)")" && pwd):/workdir" \
		-w /workdir \
		$(IMAGE_NAME) \
		python3 /app/Phylo-MIP.py "/workdir/$$(basename "$(INPUT)")" $(ARGS)
	@echo "✅ Phylo-MIP analysis complete!"

# Merge data command
merge-data:
	@if [ -z "$(QIIME)" ] || [ -z "$(PHYLOMIP)" ] || [ -z "$(FORMAT)" ]; then \
		echo "❌ Error: QIIME, PHYLOMIP, and FORMAT parameters are required"; \
		echo "Usage: make merge-data QIIME=<file> PHYLOMIP=<file> FORMAT=<format> [OUTPUT=<file>] [DEBUG=1]"; \
		echo "Example: make merge-data QIIME=qiime.tsv PHYLOMIP=phylo.csv FORMAT=csv"; \
		exit 1; \
	fi
	@if [ ! -f "$(QIIME)" ]; then \
		echo "❌ Error: QIIME file '$(QIIME)' not found"; \
		exit 1; \
	fi
	@if [ ! -f "$(PHYLOMIP)" ]; then \
		echo "❌ Error: Phylo-MIP file '$(PHYLOMIP)' not found"; \
		exit 1; \
	fi
	@echo "🔄 Running merge_data..."
	@echo "📊 QIIME file: $(QIIME)"
	@echo "🧬 Phylo-MIP file: $(PHYLOMIP)"
	@echo "📋 Format: $(FORMAT)"
	$(if $(OUTPUT),@echo "📝 Output file: $(OUTPUT)")
	$(if $(DEBUG),@echo "🐛 Debug mode: enabled")
	@docker run --rm -i $(PLATFORM_FLAG) $(USER_FLAG) \
		-v "$$(pwd):/workdir" \
		-v "$$(cd "$$(dirname "$(QIIME)")" && pwd):/input1:ro" \
		-v "$$(cd "$$(dirname "$(PHYLOMIP)")" && pwd):/input2:ro" \
		-w /workdir \
		$(IMAGE_NAME) \
		bash -c " \
			cp /input1/$$(basename "$(QIIME)") /workdir/qiime_temp && \
			cp /input2/$$(basename "$(PHYLOMIP)") /workdir/phylomip_temp && \
			python3 /app/merge_data.py -q /workdir/qiime_temp -m /workdir/phylomip_temp -f $(FORMAT) \
				$(if $(OUTPUT),-o $(OUTPUT)) \
				$(if $(DEBUG),-d) && \
			rm -f /workdir/qiime_temp /workdir/phylomip_temp \
		"
	@echo "✅ Data merging complete!"

# Check if Docker is available
check-docker:
	@if ! command -v docker >/dev/null 2>&1; then \
		echo "❌ Error: Docker is not installed or not in PATH"; \
		echo "Please install Docker and try again."; \
		echo "Visit: https://docs.docker.com/get-docker/"; \
		exit 1; \
	fi