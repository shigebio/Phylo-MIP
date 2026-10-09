#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2024-2026 <actual copyright holder(s)>
# SPDX-License-Identifier: GPL-3.0-only
#
# This file is part of Phylo-MIP.
# See the LICENSE file in the project root for the full license text.

# リポジトリ内launcherをcanonicalとしてDocker環境を構築するsetup script
# Build the Docker environment while keeping repository launchers canonical.

REPOSITORY_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

for launcher in phylo-mip merge_data; do
    if [[ ! -f "$REPOSITORY_DIR/$launcher" ]]; then
        echo "Error: repository launcher not found: $REPOSITORY_DIR/$launcher"
        exit 1
    fi
done

warn_legacy_launcher() {
    local launcher_name="$1"
    local legacy_path="${HOME:-}/bin/$launcher_name"

    if [[ -n "${HOME:-}" ]] && { [[ -e "$legacy_path" ]] || [[ -L "$legacy_path" ]]; }; then
        echo "WARNING: legacy launcher detected: $legacy_path"
        echo "WARNING: it may take precedence over the repository launcher."
        echo "WARNING: setup.sh will not delete or overwrite this file."
        echo "Inspect with: type -a $launcher_name"
    fi
}

if [[ -n "${HOME:-}" ]]; then
    warn_legacy_launcher "phylo-mip"
    warn_legacy_launcher "merge_data"
else
    echo "WARNING: HOME is not set; legacy launcher detection was skipped."
fi

# OS判定 / Detect the supported host operating system.
case "$(uname -s)" in
    Darwin*)
        OS_TYPE="macos"
        ;;
    Linux*)
        OS_TYPE="linux"

        # Linux launcher実行に必要なrealpathを確認 / Check realpath required by the Linux launcher.
        if ! command -v realpath >/dev/null 2>&1; then
            echo "Installing realpath..."
            sudo apt-get update && sudo apt-get install -y coreutils || {
                echo "Error: Failed to install coreutils."
                exit 1
            }
        fi
        ;;
    *)
        echo "Unsupported OS: $(uname -s)"
        exit 1
        ;;
esac

echo "Setting up Phylo-MIP for ${OS_TYPE}..."
echo "Repository launchers:"
echo "  $REPOSITORY_DIR/phylo-mip"
echo "  $REPOSITORY_DIR/merge_data"

for launcher in phylo-mip merge_data; do
    if [[ ! -x "$REPOSITORY_DIR/$launcher" ]]; then
        echo "WARNING: $REPOSITORY_DIR/$launcher is not executable."
        echo "Run: chmod +x ./phylo-mip ./merge_data"
    fi
done

# Keep setup side-effect free with respect to shell configuration.
# shell設定（.bashrc/.zshrc等）と$HOME/binは変更しない / Do not modify shell rc files or $HOME/bin.

# リポジトリをbuild contextに固定 / Use the repository as the Docker build context.
cd "$REPOSITORY_DIR" || exit 1
echo "Building Docker image..."
if [[ "$OS_TYPE" == "macos" ]]; then
    # Apple Siliconではlinux/amd64を指定 / Use linux/amd64 on Apple Silicon.
    if [[ "$(uname -m)" == "arm64" ]]; then
        echo "Building for ARM64 MacOS (linux/amd64 image)..."
        docker build --platform linux/amd64 -t phylo-mip "$REPOSITORY_DIR"
    else
        echo "Building for Intel MacOS..."
        docker build -t phylo-mip "$REPOSITORY_DIR"
    fi
else
    docker build -t phylo-mip "$REPOSITORY_DIR"
fi

echo "Setup complete."
echo "Canonical commands from the repository root:"
echo "  ./phylo-mip <input_file> [options]"
echo "  ./merge_data -q <qiime_file> -p <Phylo-MIP_file> -f <format> [options]"
echo "No launcher copies or shell PATH entries were created or updated."
