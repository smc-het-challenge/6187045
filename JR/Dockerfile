FROM ubuntu

RUN apt-get update; \
    apt-get install -y \
      python python-pip \
      python-numpy python-scipy \
      build-essential python-dev python-setuptools \
      libatlas-dev libatlas3gf-base

RUN update-alternatives --set libblas.so.3 \
      /usr/lib/atlas-base/atlas/libblas.so.3; \
    update-alternatives --set liblapack.so.3 \
      /usr/lib/atlas-base/atlas/liblapack.so.3

RUN pip install -U scikit-learn

# install python
#RUN apt-get update
#RUN apt-get install -y python python-dev python-distribute python-pip

# install python packages contain in requirements
RUN mkdir -p /usr/src/app
WORKDIR /usr/src/app
COPY requirements.txt /usr/src/app/
RUN pip install -r requirements.txt
COPY . /usr/src/app

# add our own scripts 
WORKDIR /opt/JR
ADD MOB.py /opt/JR/MOB.py
ADD merge_vcfBat.py /opt/JR/merge_vcfBat.py
