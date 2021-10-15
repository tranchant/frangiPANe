ARG BASE_CONTAINER=jupyter/datascience-notebook
FROM $BASE_CONTAINER

LABEL author="frangiPANe"

USER root

RUN apt update \
&&  apt install -y unzip wget build-essential cmake git-all tar gzip rsync

ENV JUPYTER_ENABLE_LAB=yes

RUN mkdir /data/ \
&& git clone  --recursive https://github.com/tranchant/frangiPANe.git /data/ \
&& ln -s /data/frangiPANe/ .

USER $NB_UID
