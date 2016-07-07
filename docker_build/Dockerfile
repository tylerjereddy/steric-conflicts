FROM ubuntu:14.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
	git \
	libfontconfig1 \
	libglib2.0-0 \
	libsm6 \
	libxrender1 \
	python-pip \
	wget 

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh && \
    /bin/bash /Miniconda2-4.0.5-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda2-4.0.5-Linux-x86_64.sh

ENV PATH /opt/conda/bin:$PATH

RUN conda install numpy=1.11.0 scipy=0.17.0 cython nose matplotlib
RUN git clone https://github.com/MDAnalysis/mdanalysis.git /src/mdanalysis
RUN pip install -e /src/mdanalysis/package
RUN pip install coveralls

ADD *py /steric_analysis/
ADD *gro /steric_analysis/

#coveralls currently requires a git repo & unique commit to submit the report, so
#do this artificially:
RUN cd /steric_analysis && \
    git config --global user.email "dummy@dummy.com" && \
    git config --global user.name "Dummy" && \
    git init && \
    git add * && \
    git commit -m "Dummy."
