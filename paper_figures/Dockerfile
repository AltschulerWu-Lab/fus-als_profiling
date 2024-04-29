## Base docker image 
FROM ubuntu:latest 
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && \
    apt-get install -y r-base \
      libssl-dev \
      libcurl4-gnutls-dev \
      libxml2-dev \
      libfontconfig1-dev \
      libharfbuzz-dev \
      libfribidi-dev \
      libfreetype6-dev \
      libpng-dev \
      libtiff5-dev \
      libjpeg-dev \
      libgit2-dev \
      git \
      curl \
      python3-venv

RUN apt install -y cmake
## TODO: install python for keras

## Initialize Renviron varaibles
RUN touch ~/.Renviron
RUN echo "ALS_PAPER=/ALS/" >> ~/.Renviron
RUN echo "ALS_DATA=/ALS/data_profiles/" >> ~/.Renviron

## run the script
COPY ./setup.R /setup.R
RUN Rscript /setup.R

