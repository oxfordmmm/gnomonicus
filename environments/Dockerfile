FROM ubuntu:20.04

#Install requirements
RUN apt update && apt install -y python3-pip git

#Install gnomon from git
RUN git clone https://github.com/oxfordmmm/gnomon && cd gnomon && pip install .
