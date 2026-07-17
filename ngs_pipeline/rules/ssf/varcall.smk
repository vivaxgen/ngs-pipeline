
from ngs_pipeline.rules import pkg

include: pkg("ngs_pipeline::ssf/init.smk")

include: pkg("ngs_pipeline::workflows/discovery-varcall.smk")