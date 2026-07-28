# Third-party software notices

Phylo-MIP installs or invokes third-party software.

These components are not relicensed under the Phylo-MIP license.
Each component remains subject to its upstream license.

| Component | Purpose | Upstream license | Source |
|---|---|---|---|
| MAFFT | Multiple-sequence alignment | BSD license | https://mafft.cbrc.jp/alignment/software/ |
| VSEARCH | Haplotype clustering | GPL-3.0-or-later OR BSD-2-Clause | https://github.com/torognes/vsearch |
| FastTree | Phylogenetic tree construction | GPL-2.0-or-later; verify the license of the installed Ubuntu package | https://github.com/morgannprice/fasttree |
| PTP / bPTP | Species-delimitation analysis | GNU GPL version 3; verify whether the installed revision is `only` or `or-later` | https://github.com/zhangjiajie/PTP |
| mPTP | Species-delimitation analysis | AGPL-3.0-or-later | https://github.com/Pas-Kapli/mptp |

Python packages installed through `requirements.txt` and operating-system
packages installed through `apt` are distributed under their own licenses.

The complete list of installed packages and versions depends on the versions
resolved when the Docker image is built.

Before distributing a pre-built Docker image, the distributor must verify the
licenses of the exact installed versions, preserve required copyright and
license notices, and satisfy any corresponding-source requirements.