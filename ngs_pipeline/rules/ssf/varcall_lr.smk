# SPDX-FileCopyrightText: 2026 Hidayat Trimarsanto <trimarsanto@gmail.com>
# SPDX-License-Identifier: MIT

__copyright__ = "(c) 2026 Hidayat (Anto) Trimarsanto <trimarsanto@gmail.com>"
__author__ = "trimarsanto@gmail.com"
__license__ = "MIT"

from ngs_pipeline.rules import pkg

# set the reads file name as being set by ngs-pl prepare-sample-directory command
config["reads_file"] = config.get("reads_file", "reads/raw-{idx}_R0.fastq.gz")

# use null trimmer since we rely on mappers to perform soft-clipping
# on primers and adapters:
config["trimmer_wf"] = config.get("trimmer_wf", "ngs_pipeline::trimmer/fastplong.smk")

# for panel sequencing, we default to minimap2
config["mapper_wf"] = config.get("mapper_wf", "ngs_pipeline::mapper/minimap2_lr.smk")

# for panel variant calling, we default to freebayes
config["varcaller_wf"] = config.get("varcaller_wf", "ngs_pipeline::varcaller/clair3.smk")

include: pkg("ngs_pipeline::ssf/varcall.smk")


# EOF
