FROM python:3.10

#Install requirements
RUN apt update && apt install -y git

#Install gnomonicus from git (for now)
# RUN git clone https://github.com/oxfordmmm/gnomonicus && cd gnomonicus && pip install .

#Install gnomonicus from pypi
RUN pip install --upgrade gnomonicus --break-system-packages
