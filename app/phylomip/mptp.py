"""mPTP execution and mPTP output directory handling."""

import os
import subprocess


def run_mptp(tree_file, base_dir, timestamp, subprocess_module=subprocess):
    try:
        output_dir = os.path.join(base_dir, f"{timestamp}_mPTP_analysis")
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f"{timestamp}_mPTP_species_delimitation")
        command = ["xvfb-run", "-a", "mptp", "-tree_file", tree_file,
                   "-output_file", output_file, "-ml", "-single"]
        subprocess_module.run(command, check=True)
        print("mPTP analysis complete.")
    except subprocess_module.CalledProcessError as exc:
        print(f"Error running mPTP: {exc}")

