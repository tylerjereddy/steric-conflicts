FROM ubuntu:14.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
	git \
	python-pip \
	wget 

RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh && \
    /bin/bash /Miniconda2-4.0.5-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda2-4.0.5-Linux-x86_64.sh

ENV PATH /opt/conda/bin:$PATH
#shouldn't expose this publicly, but for now seems necessary to report coveralls
#result from INSIDE the docker container INSIDE Travis-CI:
ENV COVERALLS_REPO_TOKEN tqiadfk3Gzqk4HDdoozi71Fg7irjQBKkD

RUN conda install numpy scipy
RUN conda config --add channels MDAnalysis && conda install mdanalysis=0.14
RUN pip install coveralls

ADD *py /steric_analysis/
ADD *gro /steric_analysis/

#coveralls currently requires a git repo & unique commit to submit the report, so
#do this artificially:
RUN cd /steric_analysis && git config --global user.email "dummy@dummy.com" && \
     git config --global user.name "Dummy" && git init && git commit -m "Testing."
