ARG BASE_CONTAINER=jupyter/datascience-notebook
FROM $BASE_CONTAINER

LABEL author="frangiPANe"

USER root

RUN apt-get update \
&&  apt-get install -y less libz-dev software-properties-common apt-utils unzip wget build-essential cmake git-all tar gzip rsync

RUN pip install bash_kernel
RUN python3 -m bash_kernel.install

RUN pip3 install panel \
&& pip3 install biopython

RUN curl -fsSL https://deb.nodesource.com/setup_16.x | sudo -E bash - \
&& sudo apt-get install -y nodejs

RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager

ENV JUPYTER_ENABLE_LAB=yes

RUN mkdir /mydata/ && chown $NB_USER /mydata -R

RUN mkdir /data/ \
&& git clone  --recursive https://github.com/tranchant/frangiPANe.git /data/ \
&& ln -s /data/frangiPANe/ /home/jovyan/

RUN apt-get install -y ea-utils bwa samtools
RUN apt-get install -y abyss cd-hit
RUN apt-get install -y ncbi-blast+ ncbi-tools-bin

RUN conda install -c bioconda assembly-stats

USER $NB_UID
#  https://stackoverflow.com/questions/55362117/symlink-in-docker-container-not-supported


#RUN mkdir /home/jovyan/mydata

# Create a mountpoint
#RUN mkdir /mydata
