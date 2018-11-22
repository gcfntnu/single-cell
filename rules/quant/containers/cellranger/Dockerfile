FROM ubuntu:18.04

ENV CELLRANGER_VER 3.0.0
ENV BCL2FASTQ_VER 2-20-0

RUN apt-get update --fix-missing && \
    apt-get install -y \
    alien \
    unzip \
    wget 

RUN \ 
  wget http://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v$BCL2FASTQ_VER-linux-x86-64.zip && \
  unzip bcl2fastq2-v$BCL2FASTQ_VER-linux-x86-64.zip  && \
  alien bcl2fastq2*.rpm && \
  dpkg -i bcl2fastq2*.deb && \
  rm bcl2fastq2*.deb bcl2fastq2*.rpm bcl2fastq2*.zip

COPY files/cellranger-$CELLRANGER_VER.tar.gz /opt/cellranger-$CELLRANGER_VER.tar.gz
RUN \
  cd /opt && \
  tar -xzvf cellranger-$CELLRANGER_VER.tar.gz && \
  export PATH=/opt/cellranger-$CELLRANGER_VER:$PATH && \
  ln -s /opt/cellranger-$CELLRANGER_VER/cellranger /usr/bin/cellranger && \
  rm -rf /opt/cellranger-$CELLRANGER_VER.tar.gz
