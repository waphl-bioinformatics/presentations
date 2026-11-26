#!/bin/bash

rm -r results bigbacter_db .nextflow* work/
nextflow run ../../bigbacter/main.nf -profile docker,pc -params-file params/params.1.yml
nextflow run ../../bigbacter/main.nf -profile docker,pc -params-file params/params.2.yml
nextflow run ../../bigbacter/main.nf -profile docker,pc -params-file params/params.3.yml