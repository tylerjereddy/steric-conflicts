FROM ubuntu:14.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
	git \
	wget 

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh && \
    /bin/bash /Miniconda2-4.0.5-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda2-4.0.5-Linux-x86_64.sh

ENV PATH /opt/conda/bin:$PATH

RUN conda install numpy scipy
RUN conda config --add channels MDAnalysis && conda install mdanalysis=0.14

ADD *py /steric_analysis/
ADD *gro /steric_analysis/
