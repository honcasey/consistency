# hash:sha256:ed3af7e148046c1a06a8aaf51c19c9eb62179400c8d2771d53e7a06dcd103735
FROM registry.codeocean.com/codeocean/r-base:4.0.0-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys \
        0xAD2323F17326AE31401037733E05EBFF05441C52 \
    && add-apt-repository -y 'deb http://deb.codeocean.com/rstudio-server-bionic/ ubuntu main' \
    && apt-get purge -y --autoremove software-properties-common \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        pandoc=1.19.2.4~dfsg-1build4 \
        rstudio-server=1.2.5033 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("VennDiagram", "1.6.20")' \
    && Rscript -e 'remotes::install_version("dplyr", "1.0.1")' \
    && Rscript -e 'remotes::install_version("survival", "3.2-3")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "PharmacoGx", \
        "marray" \
    ))' # Original versions: 2.0.5 1.66.0
