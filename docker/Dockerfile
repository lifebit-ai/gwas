FROM lifebitai/bcftools:1.9
LABEL authors="phil@lifebit.ai" \
      description="Docker image containing dependencies for GEL GWAS"

RUN apt-get update && \
    apt-get install unzip -y
RUN wget http://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip && \
    unzip plink2_linux_x86_64.zip -d plink2 && \
    rm plink2_linux_x86_64.zip

ENV PATH="$PATH:/plink2"