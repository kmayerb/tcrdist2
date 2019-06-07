FROM ubuntu:18.04

RUN apt-get update
RUN apt-get -y install git python python-pip libpng-dev libxft-dev imagemagick wget

# install python dependencies
RUN pip install numpy==1.10.1
RUN pip install scipy==0.16.0
RUN pip install scikit-learn==0.17
RUN pip install matplotlib==1.4.3

ADD Legacy /tcrdist2

WORKDIR /tcrdist2
# set the environment variable that we are working on a docker build
ENV DOCKER_BUILD 1
RUN python setup.py 