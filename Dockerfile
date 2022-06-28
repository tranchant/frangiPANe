ARG BASE_CONTAINER=jupyter/datascience-notebook
FROM $BASE_CONTAINER

LABEL author="frangiPANe"

USER root

RUN apt-get update \
&&  apt-get install -y less libz-dev software-properties-common apt-utils unzip wget build-essential cmake git-all tar gzip rsync python3-pip

RUN apt pip3 install --upgrade pip
# RUN sudo pip3 install bash_kernel
# RUN sudo python3 -m bash_kernel.install

RUN pip3 install panel \
&& pip3 install biopython

RUN curl -fsSL https://deb.nodesource.com/setup_16.x | sudo -E bash - \
&& apt-get install -y nodejs

RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager

ENV JUPYTER_ENABLE_LAB=yes

RUN mkdir /data/ \
&& git clone  --recursive https://github.com/tranchant/frangiPANe.git /data/ \
&& ln -s /data/frangiPANe/ /home/jovyan/ && sudo chown $NB_USER /data/frangiPANe/ -R

RUN apt-get install -y ea-utils bwa samtools
RUN apt-get install -y abyss cd-hit
RUN apt-get install -y ncbi-blast+ ncbi-tools-bin

RUN conda install -c bioconda assembly-stats

RUN mkdir /mydata/ && sudo chown $NB_USER /mydata -R

USER $NB_UID
#  https://stackoverflow.com/questions/55362117/symlink-in-docker-container-not-supported


#RUN mkdir /home/jovyan/mydata

# Create a mountpoint
#RUN mkdir /mydata
