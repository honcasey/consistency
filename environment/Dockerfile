# hash:sha256:b0d01c356d7e1fe766b20d1e6f20554de0809cf7bc252dab67bc7355a716719c
FROM registry.codeocean.com/codeocean/r-base:3.6.0-ubuntu18.04

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
    && Rscript -e 'remotes::install_version("survival", "3.2-3")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "PharmacoGx", \
        "marray" \
    ))' # Original versions: 1.17.1 1.64.0
