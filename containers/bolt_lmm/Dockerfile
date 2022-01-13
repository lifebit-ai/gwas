#quay.io/lifebitai/bolt-lmm:2.3.6

FROM continuumio/miniconda3@sha256:456e3196bf3ffb13fee7c9216db4b18b5e6f4d37090b31df3e0309926e98cfe2
LABEL authors="eva@lifebit.ai" \
      description="Docker image containing BOLT-LMM"


RUN apt-get update --allow-releaseinfo-change && \
    apt-get install -y \
    procps && \
    apt-get clean -y 

RUN wget https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.3.6.tar.gz \
  && tar -xvf BOLT-LMM_v2.3.6.tar.gz \
  && rm BOLT-LMM_v2.3.6.tar.gz

ENV PATH="$PATH:/BOLT-LMM_v2.3.6/"

ENTRYPOINT []