FROM ubuntu:latest

#Install requirements
RUN apt update && apt install -y git python3 python3-pip

#Install gnomonicus from git
# RUN git clone https://github.com/oxfordmmm/gnomonicus && cd gnomonicus && pip install .

#Install gnomonicus from pypi
RUN pip install --upgrade gnomonicus
