FROM ubuntu
MAINTAINER joseph.schaeffer@autodesk.com

# USAGE: To build this instance:
#        docker build -t <identifier> .
#      where <identifier> is what you use as an identifier for the resulting container. Typically, 
#      a username / project name combination is reasonable for use here. 
# 
#        To run this instance (and be placed in a bash shell):
#        docker run -t -i <identifier> /bin/bash
#


RUN apt-get update -y --fix-missing
RUN apt-get upgrade -y
RUN apt-get install python2.7 python2.7-dev python-pip -y
# May need git or wget here for future dependencies.

# Any pip installs can go here. Currently no requirements.
# RUN pip install --upgrade biopython
RUN pip install --upgrade numpy

ENV MODULE /nanodesign

RUN mkdir -p $MODULE

COPY ./ ${MODULE}/

# install the module
RUN cd $MODULE && python setup.py install

# set up the scripts directory
ENV APP /app
WORKDIR $APP

COPY scripts/ ${APP}/
