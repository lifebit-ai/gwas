#quay.io/lifebitai/regenie

FROM ubuntu:18.04
LABEL authors="eva@lifebit.ai" \
      description="Docker image containing Regenie"


RUN apt-get update --allow-releaseinfo-change && \
    apt-get install -y \
    wget \
    libgomp1 \
    zip \
    procps && \
    apt-get clean -y 

RUN wget https://github.com/rgcgithub/regenie/releases/download/v2.2.4/regenie_v2.2.4.gz_x86_64_Linux.zip

RUN mkdir -p REGENIE/

RUN unzip regenie_v2.2.4.gz_x86_64_Linux.zip -d REGENIE/

RUN mv REGENIE/regenie_v2.2.4.gz_x86_64_Linux REGENIE/regenie

ENV PATH="$PATH:/REGENIE/"

ENTRYPOINT []