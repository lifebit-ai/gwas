FROM nfcore/base:1.10.2
LABEL authors="Marcos Camara Donoso" \
      description="Docker image containing all software requirements for the ldsc"

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

RUN git clone https://github.com/bulik/ldsc.git && \
    cd ldsc && \
    conda env create --quiet -f environment.yml && conda clean -a

# Correct the shebang of the files to point to python2
RUN for i in `ls /ldsc/*py` ; do sed -i 's/python/python2/g'  $i; done
RUN for i in `ls /ldsc/ldscore/*py | grep -v __init__` ; do sed -i '1 i #!/usr/bin/env python2' $i ; done

ENV PATH /opt/conda/envs/ldsc/bin:$PATH

RUN mv /opt/conda/envs/ldsc/bin/python /opt/conda/envs/ldsc/bin/python3
RUN mv /opt/conda/envs/ldsc/bin/python2 /opt/conda/envs/ldsc/bin/python

RUN find /ldsc/ -type f -iname "*.py" -exec chmod x {} \;
RUN touch .Rprofile
RUN touch .Renviron

#Add assets/ folder
RUN mkdir /assets/
COPY assets/* /assets/

ENV PATH="$PATH:/ldsc/"