#!/bin/bash
##############################################################################################################################
# pipeline_template
# 
# This program is free software that contains third party software subject to various licenses, 
# namely, the GNU General Public License version 3 (or later), the GNU Affero General Public License 
# version 3 (or later), and the LaTeX Project Public License v.1.3(c). A list of the software contained 
# in this program, including the applicable licenses, can be accessed here: 
# 
#
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
#
# Program:  download-genome-fa.sh
# Version:  v0.1
# Author:   Sami R. Cherikh
# Purpose:  Download FASTA genome from NCBI ftp
# Input:    N/A
# Output:   N/A
##############################################################################################################################

## Command line arguments
resDir=$1

## migrate to results directory
cd $resDir;

## get Genome fa from NCBI ftp
echo "downloading Genome to map reads to GRCh38 or hg38..."
echo "wget \"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz\""
wget -nv "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz";

## Unzip
gunzip -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz;

echo 'Done with genome...'