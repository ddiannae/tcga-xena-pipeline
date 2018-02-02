FROM r-base

MAINTAINER Diana Garcia <diana.gco@gmail.com>

ADD . /pipeline
WORKDIR /pipeline
ENV HOME /pipeline

RUN apt-get update -qq && apt-get install -y \
	wget \
    unzip

RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("edgeR");'

#Install GDC Data Transfer Tool Client
RUN wget http://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.3.0_Ubuntu14.04_x64.zip
RUN unzip gdc-client_v1.3.0_Ubuntu14.04_x64.zip
RUN chmod 755 gdc-client && ln -s /pipeline/gdc-client /usr/local/bin/

RUN rm *.zip

ADD start.sh /
CMD ["/start.sh"]

