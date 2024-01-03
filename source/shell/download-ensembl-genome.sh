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
# Program:  download-ensembl-genome.sh
# Version:  v0.1
# Author:   Travis L. Jensen and Johannes B. Goll
# Purpose:  Download FASTA formatted Ensembl genome.  Please note chromosomes are NOT ordered
#				in the resulting genome.
# Input:    N/A
# Output:   N/A
##############################################################################################################################

## Command line arguments
## $1) ensembl version # ie: 87
## $2) download results to directory ie: /home/user01/annot
ensemblVer=$1
resDir=$2

## migrate to results directory
cd $resDir;

## get Ensembl Genome
echo "downloading Ensembl Version $ensemblVer Genome..."
echo "wget \"ftp://ftp.ensembl.org/pub/release-$ensemblVer/fasta/homo_sapiens/dna/Homo_sapiens.*.dna.chromosome.*.fa.gz\""
wget -nv "ftp://ftp.ensembl.org/pub/release-$ensemblVer/fasta/homo_sapiens/dna/Homo_sapiens.*.dna.chromosome.*.fa.gz";

## merge chromosomes
echo "zcat Homo_sapiens.*.dna.chromosome.*.fa.gz > Homo_sapiens.ensembl.version$ensemblVer.genome.fa";
zcat Homo_sapiens.*.dna.chromosome.*.fa.gz > Homo_sapiens.ensembl.version$ensemblVer.genome.fa;

# clean up
rm *.fa.gz;

echo 'Done with genome...'