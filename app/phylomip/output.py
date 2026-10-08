"""Output directory layout, naming, copying, and verification."""

import os
import shutil
import stat
from dataclasses import dataclass
from datetime import datetime


@dataclass(frozen=True)
class OutputPaths:
    input_dir: str
    output_dir: str
    timestamp: str
    taxonomy_dir: str
    alignment_dir: str
    phylogeny_dir: str
    bptp_base_dir: str
    mptp_base_dir: str

    @classmethod
    def from_input(cls, input_csv, timestamp=None):
        input_file = os.path.abspath(input_csv)
        input_dir = os.path.dirname(input_file)
        stamp = timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = os.path.join(input_dir, f"phylomip_output_{stamp}")
        return cls(input_dir, output_dir, stamp,
                   os.path.join(output_dir, "taxonomy"),
                   os.path.join(output_dir, "alignment"),
                   os.path.join(output_dir, "phylogeny"),
                   os.path.join(output_dir, "bptp"),
                   os.path.join(output_dir, "mptp"))


def ensure_directory_exists(file_path):
    directory = os.path.dirname(str(file_path))
    if not directory or os.path.exists(directory):
        return
    try:
        os.makedirs(directory, exist_ok=True)
        try:
            os.chmod(directory, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
        except (OSError, PermissionError):
            pass
    except Exception as exc:
        print(f"Warning: Could not create directory {directory}: {exc}")


def create_output_directories(paths):
    for directory in (paths.output_dir, paths.taxonomy_dir, paths.alignment_dir,
                      paths.phylogeny_dir, paths.bptp_base_dir, paths.mptp_base_dir):
        try:
            os.makedirs(directory, exist_ok=True)
            try:
                os.chmod(directory, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
            except (OSError, PermissionError) as exc:
                print(f"Warning: Could not set permissions for {directory}: {exc}")
            print(f"Created directory: {directory}")
            if not os.path.exists(directory):
                print(f"ERROR: Failed to create directory: {directory}")
            else:
                print(f"SUCCESS: Directory exists: {directory}")
        except Exception as exc:
            print(f"Error creating directory {directory}: {exc}")
            raise


def verify_directory_structure(paths):
    print("\n=== Directory Structure Verification ===")
    for root, dirs, files in os.walk(paths.output_dir):
        level = root.replace(paths.output_dir, "").count(os.sep)
        indent = "  " * level
        print(f"{indent}{os.path.basename(root)}/")
        for file_name in files:
            print(f"{'  ' * (level + 1)}{file_name}")
    print("=========================================\n")


def copy_with_directory_structure(src_dir, dest_dir):
    try:
        if os.path.exists(dest_dir):
            shutil.rmtree(dest_dir)
        shutil.copytree(src_dir, dest_dir)
        print(f"Successfully copied directory structure from {src_dir} to {dest_dir}")
    except Exception as exc:
        print(f"Error copying directory structure: {exc}")


def copy_output_files_with_structure(output_dir, original_dir):
    output_dirname = os.path.basename(output_dir)
    dest_path = os.path.join(original_dir, f"copied_{output_dirname}")
    try:
        if os.path.exists(dest_path):
            shutil.rmtree(dest_path)
        shutil.copytree(output_dir, dest_path)
        print(f"Successfully copied directory to: {dest_path}")
        return dest_path
    except Exception as exc:
        print(f"Error copying output files: {exc}")
        return None

