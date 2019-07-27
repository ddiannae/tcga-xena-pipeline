FROM rocker/r-ver:3.4.4 AS builder-r
RUN apt-get update && \
  apt-get -y install --fix-missing libssl-dev libcurl4-openssl-dev libxml2-dev zlib1g-dev libmariadb-client-lgpl-dev libpng-dev 
COPY 00_InstallPackages.R . 
RUN R -f 00_InstallPackages.R

FROM rocker/r-ver:3.4.4
WORKDIR /pipeline
COPY --from=builder-r /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY ["Biomart_EnsemblG94_GRCh38_p12.txt", "01_GetTheData.R", "02_QC.R", "03_Normalization.R", "start.sh", "/pipeline/"]
COPY gdc-client /usr/local/bin/
RUN apt-get update \
  && apt-get -y install --fix-missing libssl-dev libcurl4-openssl-dev libxml2-dev zlib1g-dev libmariadb-client-lgpl-dev libpng-dev 
RUN useradd -m pipeuser -u 1033 \
  &&  groupmod -g 1003 $(id pipeuser -g -n) \
  && chown -R pipeuser /pipeline
USER pipeuser
CMD ["/pipeline/start.sh"]
