#!/bin/bash

# to be sourced
DB=meta4_bins
HOST=`mysql_host`

# pathes
BWA=/home/dacb/express/software/bwa/bin/bwa
SAMTOOLS=/home/dacb/express/software/samtools/bin/samtools
MYSQL=mysql
MYSQL_HOST=/home/dacb/bin/mysql_host
HTSEQ_COUNT=/home/dacb/express/software/htseq/bin/htseq-count

# job parameters
QL=walltime=999:99:99,mem=2gb,feature=8core
GROUP_LIST=hyak-lidstrom
EMAIL=dacb@uw.edu,jmatsen@uw.edu
