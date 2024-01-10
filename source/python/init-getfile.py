#############################################################################################################
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
# Program:  init-getfile.py
# Version:  0.1
# Author:   Sami R. Cherikh, Travis L. Jensen
# Purpose:  Initializer script for pulling input files
# Input:    N/A
# Output:   N/A
#############################################################################################################

## Import modules
import sys, argparse

parser=argparse.ArgumentParser()

parser.add_argument("--src", help="Source dir")
parser.add_argument("--ends", help="Single or paired end indication")
parser.add_argument("--samid", help="Sample IDs")
parser.add_argument("--fastq1", help="Fastq read file 1", default='', required=False, nargs='?', const='')
parser.add_argument("--fastq2", help="Fastq read file 2", default='', required=False, nargs='?', const='')
parser.add_argument("--bam", help="Bam read input", default='', required=False, nargs='?', const='')
parser.add_argument("--cloud", help="Cloud program name")
parser.add_argument("--fastqdump", help="SRA kit fastqdump program name")
parser.add_argument("--openssl", help="Openssl program name", default='', required=False, nargs='?', const='')
parser.add_argument("--decrypt-pass", help="Decryption pswd", default='', required=False, nargs='?', const='')
parser.add_argument("--hash", help="Hash alg", default='', required=False, nargs='?', const='')
parser.add_argument("--sample", help="Sample identifier")
parser.add_argument("--output", help="Output")

args=parser.parse_args()


## Use the source dir to import helper modules
sys.path.append(args.src+'/python')
import getfile
import utils


## process args
ENDS = str(args.ends).split(',')
SAMID = str(args.samid).split(',')
FASTQ_1 = str(args.fastq1).split(',')
FASTQ_2 = str(args.fastq2).split(',')
BAM = str(args.bam).split(',')
CLOUD = args.cloud
FASTQDUMP = args.fastqdump
OPENSSL = args.openssl
DECRYPT_PASS = args.decrypt_pass
HASH = args.hash

sample = args.sample
output = str(args.output).split(',')


## Get files to pull
if ENDS == ['1']:
    in_file = [FASTQ_1[i[0]] for i in enumerate(SAMID) if SAMID[i[0]] == sample]
elif ENDS != ['1'] and BAM != [''] * len(SAMID):
    in_file = [BAM[i[0]] for i in enumerate(SAMID) if SAMID[i[0]] == sample]
else:
    in_file_1 = [FASTQ_1[i[0]] for i in enumerate(SAMID) if SAMID[i[0]] == sample][0]
    in_file_2 = [FASTQ_2[i[0]] for i in enumerate(SAMID) if SAMID[i[0]] == sample][0]
    in_file = [in_file_1, in_file_2]
    ## if using SRA that are paired end, only pass first file but will use split-3 in fastq-dump call
    if(in_file_1[0][0:3] == 'SRR'):
        in_file = [in_file_1]


## Run command
for f in range(len(in_file)):
    ## Call getFile for each in file
    getfile.getFile(in_file=in_file[f], out_file=utils.toList(output)[f], cloud_prog=CLOUD, fastq_dump_prog=FASTQDUMP, openssl_prog=OPENSSL, pw=DECRYPT_PASS, hash=HASH, end=len(ENDS))

