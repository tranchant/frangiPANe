ARG BASE_CONTAINER=jupyter/datascience-notebook
FROM $BASE_CONTAINER

LABEL author="frangiPANe"

USER root

RUN apt update \
&&  apt install -y libz-dev software-properties-common apt-utils unzip wget build-essential cmake git-all tar gzip rsync


RUN apt install -y python3.9-tk

# RUN jupyter labextension install @pyviz/jupyterlab_pyviz

RUN pip3 install panel \
&& pip3 install biopython


RUN curl -fsSL https://deb.nodesource.com/setup_16.x | sudo -E bash - \
&& sudo apt-get install -y nodejs

RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager


ENV JUPYTER_ENABLE_LAB=yes

RUN mkdir /data/ \
&& git clone  --recursive https://github.com/tranchant/frangiPANe.git /data/ \
&& ln -s /data/frangiPANe/ /home/jovyan/

#RUN mkdir /usr/local/bwa && cd /usr/local/bwa \
#&& git clone https://github.com/lh3/bwa.git \
#&& cd bwa; make

RUN apt install -y ea-utils bwa samtools
RUN apt install -y abyss 

ENV PATH="/usr/local/bwa/bwa:$PATH"
RUN echo $PATH

USER $NB_UID
#  https://stackoverflow.com/questions/55362117/symlink-in-docker-container-not-supported
