#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
DATA_DIR="$SCRIPT_DIR/../datasets/"

set -euxo pipefail
command -v NGmerge || conda install -y -c bioconda ngmerge
command -v flash   || conda install -y -c bioconda flash
command -v vsearch || conda install -y -c bioconda vmerge

gunzip "$DATA_DIR"/*.gz

OUT="$SCRIPT_DIR/out";
mkdir -p "$OUT"
rm -rf "$OUT"/*
set -euxo pipefail

MIN=17
MAX=299

for R1 in $DATA_DIR/*_R1*;
do
	R1=$(readlink -f "$R1")
	R2=$(readlink -f "${R1/_R1/_R2}");
	B=$(basename "$R1" | cut -f1 -d_);
	echo FOR=$R1
	echo REV=$R2
	flash   -m $MIN -M $MAX -d "$OUT" -o "$B.flash" "$R1" "$R2"
	mv "$OUT/$B.flash.extendedFrags.fastq" "$OUT/$B.flash.fq"
	NGmerge -m $MIN -o "$OUT/$B.ngmerge.fq" -1 "$R1" -2 "$R2"
	vsearch --fastq_mergepairs "$R1" --reverse "$R2" --fastqout "$OUT/$B.vsearch.fq"
done

gzip "$DATA_DIR"/*q
