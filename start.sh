#!/bin/bash
# This is a comment!
echo Hello World
echo $PRIMARY_SITE

DATADIR="data/"
install -d -m 775 data/counts
CANCERDIR="data/counts/cancer"
install -d -m 775 $CANCERDIR
HEALTHYDIR="data/counts/healthy"
install -d -m 775 $HEALTHYDIR
LOGSDIR="data/logs"
install -d -m 775 $LOGSDIR

if [[ $NODOWNLOAD == "true" ]]; then
  echo No data downloaded
else
  gdc-client download -d $CANCERDIR -m manifests/${PRIMARY_SITE}.cancer.txt --log-file ${LOGSDIR}/${PRIMARY_SITE}.cancer.log --retry-amount 3
  gdc-client download -d $HEALTHYDIR -m manifests/${PRIMARY_SITE}.healthy.txt --log-file ${LOGSDIR}/${PRIMARY_SITE}.healthy.log --retry-amount 3

  find $CANCERDIR -name '*.gz' -exec mv '{}' $CANCERDIR \;
  find $HEALTHYDIR -name '*.gz' -exec mv '{}' $HEALTHYDIR \;

  find . -type d -empty -delete
  find . -name '*.gz' -exec gunzip '{}' \;
  find . -name '*.gz' -exec rm '{}' \;

  echo Finished downloading
fi

if [[ ! $NOGETDATA == "true" ]]; then
  Rscript /pipeline/01_GetTheData.R $DATADIR
fi

if [[ ! $NOQC == "true" ]]; then
  Rscript /pipeline/02_QC.R $DATADIR
fi

if [[ ! $NONORMALIZATION == "true" ]]; then
  Rscript /pipeline/03_Normalization.R $DATADIR
fi
