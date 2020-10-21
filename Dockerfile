FROM nfcore/base:1.10.2
LABEL authors="Phil Palmer, Christina Chatzipantsiou" \
      description="Docker image containing all software requirements for the gel-gwas pipeline"

RUN apt-get update && \
    apt-get install -y \
              build-essential \
              git \
              unzip \
              autoconf \
              zlib1g-dev \
              libbz2-dev \
              liblzma-dev \
              libcurl4-gnutls-dev \
              libssl-dev \
              libgsl0-dev \
              libperl-dev \
              libxt-dev \
              speedtest-cli \
              procps

#Install htslib developmental version specific commit
RUN mkdir htslib && \
    cd htslib && \
    git init && \
    git remote add origin  git://github.com/samtools/htslib.git && \
    git fetch origin 31f0a76d338c9bf3a6893b71dd437ef5fcaaea0e && \
    git reset --hard FETCH_HEAD && \
    autoheader && \
    autoconf && \
    ./configure --enable-libcurl --enable-s3 --prefix=/ && \
    make && \
    make install && \
    cd /

#Install bcftools developmental version specific commit
RUN mkdir bcftools && \
    cd bcftools && \
    git init && \
    git remote add origin  git://github.com/samtools/bcftools.git && \
    git fetch origin 7cd83b71405c993d2f81fb91126abe3e44344398 && \
    git reset --hard FETCH_HEAD && \
    autoheader && \
    autoconf && \
    ./configure --enable-libgsl --enable-perl-filters --prefix=/ && \
    make && \
    make install && \
    cd /

RUN apt-get update && \
    apt-get install unzip -y

ENV BCFTOOLS_PLUGINS="/libexec/bcftools"

# Install plink (1)
RUN wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20200428.zip && \
    unzip plink_linux_x86_64_20200428.zip -d plink && \
    rm plink_linux_x86_64_20200428.zip

ENV PATH="$PATH:/plink"

# Install plink2
RUN wget http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip && \
    unzip plink2_linux_x86_64.zip -d plink2 && \
    rm plink2_linux_x86_64.zip

ENV PATH="$PATH:/plink2"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/gel-gwas-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name gel-gwas-1.0dev > gel-gwas-1.0dev.yml

RUN mkdir /opt/bin
COPY bin/* /opt/bin/

RUN find /opt/bin/ -type f -iname "*.py" -exec chmod +x {} \; && \
    find /opt/bin/ -type f -iname "*.R" -exec chmod +x {} \; && \
    find /opt/bin/ -type f -iname "*.sh" -exec chmod +x {} \; && \
    find /opt/bin/ -type f -iname "*.css" -exec chmod +x {} \; && \
    find /opt/bin/ -type f -iname "*.Rmd" -exec chmod +x {} \;

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

ENV PATH="$PATH:/opt/bin/"