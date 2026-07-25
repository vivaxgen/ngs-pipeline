
from ngs_pipeline.rules import pkg

include: pkg("ngs_pipeline::ssf/init.smk")

include: pkg("ngs_pipeline::workflow/discovery-varcall.smk")


# EOF
