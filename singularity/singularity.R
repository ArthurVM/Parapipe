Bootstrap: docker
From: ubuntu:20.04
%labels
maintainer="morrisa28@cardiff.ac.uk"
about.summary="parasite genomics pipeline container"

%post

export DEBIAN_FRONTEND=noninteractive

python_version=3.9.0
ujson_version=3.2.0
pymoi_version=0.1a
bcftools_version=1.12
freebayes_version=1.3.2-2

PACKAGES="cmake procps curl git build-essential wget zlib1g-dev software-properties-common libcurl4-gnutls-dev libxml2-dev libssl-dev"

apt-get update \
&& apt upgrade \
&& apt-get install -y $PACKAGES

apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
&& add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' \
&& apt update \
&& apt install -y r-base r-base-core r-recommended r-base-dev

echo ""
echo ""
echo " INSTALLING R PACKAGES "
echo ""
echo ""
Rscript -e 'install.packages(c("BiocManager", "optparse", "httr", "glue", "digest"), repos="https://cran.rstudio.com/")'
Rscript -e 'BiocManager::install(version = "3.20")'
Rscript -e 'require(BiocManager); BiocManager::install("SeqArray")'

%environment

%runscript
exec /bin/bash "$@"
%startscript
exec /bin/bash "$@"
