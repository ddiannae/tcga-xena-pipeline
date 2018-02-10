FROM r-base AS builder
WORKDIR /pipeline
COPY 00_InstallPackages.R . 
RUN R -f 00_InstallPackages.R

FROM r-base 
WORKDIR /pipeline
COPY --from=builder /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY gdc-client /usr/local/bin/
COPY ["start.sh", "Biomart_EnsemblG91_GRCh38_p10.txt", "01_GetTheData.R", "/"]

CMD ["/start.sh"]

