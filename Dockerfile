FROM ubuntu:latest

RUN apt-get update -y
RUN apt-get install build-essential -y
RUN apt-get install git -y

ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

RUN apt-get install apt-utils -y
RUN apt-get install cmake -y
RUN apt-get install mpich -y

RUN git clone https://bitbucket.org/azadcse/hipmcl.git

WORKDIR /hipmcl

RUN cmake .
RUN make
