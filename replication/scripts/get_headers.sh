#!/bin/bash

cd ../data

path_to_phesant="/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/derived/phesant_mod"
head -n 1 ${path_to_phesant}/*.txt | sed '/^$/d' > temp1
paste <(grep "^==" temp1 | sed 's/==> //g' | sed 's/ <==//g') <(grep -v "^==" temp1) > temp2
head temp2
