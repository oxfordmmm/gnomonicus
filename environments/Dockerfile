FROM ubuntu:20.04

#Install requirements
RUN apt update && apt install -y git python3 python3-pip

#Install gnomon from git
RUN git clone https://github.com/oxfordmmm/gnomon && cd gnomon && pip install .

#Install gnomon from pypi
# RUN pip install gnomon
