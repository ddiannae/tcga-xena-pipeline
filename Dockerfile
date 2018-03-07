FROM r-base AS builder
RUN apt-get update && \
apt-get -y install --fix-missing libcurl4-openssl-dev libxml2-dev
COPY 00_InstallPackages.R . 
RUN R -f 00_InstallPackages.R

FROM r-base 
WORKDIR /pipeline
COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library
RUN apt-get update && \
apt-get -y install --fix-missing libcurl4-openssl-dev libxml2-dev
COPY gdc-client /usr/local/bin/
COPY ["start.sh", "Biomart_EnsemblG91_GRCh38_p10.txt", "01_GetTheData.R", "02_QC.R", "/"]

CMD ["/start.sh"]
