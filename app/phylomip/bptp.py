"""bPTP execution and bPTP output directory handling."""

import os
import subprocess


def run_bptp(tree_file, mcmc, thinning, burnin, seed, base_dir, timestamp, subprocess_module=subprocess):
    try:
        output_dir = os.path.join(base_dir, f"{timestamp}_bPTP_analysis")
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f"{timestamp}_bPTP_species_delimitation")
        command = ["xvfb-run", "-a", "python3", "/app/PTP/bin/bPTP.py", "-t", tree_file,
                   "-o", output_file, "-s", str(seed), "-i", str(mcmc), "-n", str(thinning), "-b", str(burnin)]
        print(f"Running bPTP with command: {' '.join(command)}")
        subprocess_module.run(command, check=True)
        print("bPTP analysis complete.")
    except subprocess_module.CalledProcessError as exc:
        print(f"Error running bPTP.py: {exc}")

