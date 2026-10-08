#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2024-2026 <actual copyright holder(s)>
# SPDX-License-Identifier: GPL-3.0-only
# This file is part of Phylo-MIP.
# See the LICENSE file in the project root for the full license text.
"""Backward-compatible command-line entry point for Phylo-MIP."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from phylomip.pipeline import main


if __name__ == "__main__":
    main()
