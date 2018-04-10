#!/bin/bash
# This is a comment!
echo Hello World

echo $PRIMARY_SITE

DATADIR="data/"
CANCERDIR="${DATADIR}counts/cancer"
mkdir -p $CANCERDIR
HEALTHYDIR="${DATADIR}counts/healthy"
mkdir -p $HEALTHYDIR
LOGSDIR="data/logs"
mkdir -p $LOGSDIR

if [[ $DOWNLOAD == "true" ]]; then
  gdc-client download -d $CANCERDIR -m manifests/${PRIMARY_SITE}.cancer.txt --log-file ${LOGSDIR}/${PRIMARY_SITE}.cancer.log --retry-amount 3
  gdc-client download -d $HEALTHYDIR -m manifests/${PRIMARY_SITE}.healthy.txt --log-file ${LOGSDIR}/${PRIMARY_SITE}.healthy.log --retry-amount 3

  find $CANCERDIR -name '*.gz' -exec mv '{}' $CANCERDIR \;
  find $HEALTHYDIR -name '*.gz' -exec mv '{}' $HEALTHYDIR \;

  find . -type d -empty -delete
  find . -name '*.gz' -exec gunzip '{}' \;
  find . -name '*.gz' -exec rm '{}' \;

echo Finished downloading
fi

Rscript /01_GetTheData.R $DATADIR
Rscript /02_QC.R $DATADIR
