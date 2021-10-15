FROM jupyter/datascience-notebook
USER root

RUN apt update #\
&&  apt install -y unzip wget build-essential cmake git-all tar gzip rsync
 

#&& apt install -y python3-all-dev python3-pip python3-venv \
#&& apt install -y python3-pyqt5 pyqt5-dev-tools qttools5-dev-tools \
#&& pip install PyQt5 ete3 owlready2 pyproteinsExt ipympl jupyterlab \

#Adding dedicated kernel
#RUN python3 -m pip install --upgrade ipython
#RUN python3 -m pip install bash_kernel
#RUN python3 -m bash_kernel.install

ENV JUPYTER_ENABLE_LAB=yes

