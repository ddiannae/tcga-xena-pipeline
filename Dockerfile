FROM rocker/rstudio:3.4.4 AS builder
RUN apt-get update && \
  apt-get -y install --fix-missing libcurl4-openssl-dev libxml2-dev zlib1g-dev libmariadb-client-lgpl-dev 
COPY 00_InstallPackages.R . 
RUN R -f 00_InstallPackages.R

FROM rocker/rstudio:3.4.4
WORKDIR /pipeline
COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY ["Biomart_EnsemblG91_GRCh38_p10.txt", "01_GetTheData.R", "02_QC.R", "01_GetTheData_MinSet.R" , "03_Normalization.R", "/pipeline/"]
COPY ARSyN /pipeline/ARSyN
RUN apt-get update \
  && apt-get -y install --fix-missing libcurl4-openssl-dev libxml2-dev zlib1g-dev libmariadb-client-lgpl-dev 
