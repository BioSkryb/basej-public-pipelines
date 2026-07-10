# r_ado container
# Provides: R (ggplot2, dplyr, scales, reshape2) with the two ADO summary
#           scripts used by the germline/ADO analysis processes.
# Tag: r_ado_1.2
# Processes: SUMMARIZE_ADO_INTERVALS (summarize_ado_intervals.R),
#            CONCAT_SUMMARY_ADO_INTERVALS_LABELED (plot_summary_ado_intervals.R)
#
# NOTE: the original ado_vcfevalpost image was ubuntu:18.04 + the bionic-cran35
# apt mirror, which is now defunct — install.packages() failed silently and the
# image shipped without ggplot2. This uses rocker/r-ver (pinned Posit PPM snapshot)
# so the CRAN installs are reliable and the build fails if any package is missing.
#
# Build (x86/amd64 only):
#   docker build -t basejumper_r_ado_1.2 .

FROM rocker/r-ver:4.2.0

RUN R -e "install.packages(c('ggplot2','dplyr','scales','reshape2'))"

# Fail the build if any required package did not install.
RUN R -e "library(ggplot2); library(dplyr); library(scales); library(reshape2); cat('r_ado R stack OK\n')"

COPY ./summarize_ado_intervals.R /usr/local/bin/
COPY ./plot_summary_ado_intervals.R /usr/local/bin/

CMD ["/bin/bash"]
