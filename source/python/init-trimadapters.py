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
# Program:  init-trimadapters.py
# Version:  0.1
# Author:   Sami R. Cherikh, Travis L. Jensen
# Purpose:  Initializer script for trimming specified adapter sequences from input reads
# Input:    N/A
# Output:   N/A
#############################################################################################################

## Import modules
import sys, argparse


parser=argparse.ArgumentParser()

parser.add_argument("--src", help="Source dir")
parser.add_argument("--samid", help="Sample IDs")
parser.add_argument("--trim-fp", help="Trim 5 prime or not")
parser.add_argument("--trim-tp", help="Trim 3 prime or not")
parser.add_argument("--fp-adapters", help="5 prime adapters", default='', required=False, nargs='?', const='')
parser.add_argument("--tp-adapters", help="3 prime adapters", default='', required=False, nargs='?', const='')
parser.add_argument("--cutadapt", help="Adapter trimming software cutadapt")
parser.add_argument("--sample", help="Sample ID")
parser.add_argument("--input", help="Input")
parser.add_argument("--output", help="Output")
parser.add_argument("--threads", help="Threads")
parser.add_argument("--log", help="Log file")

args=parser.parse_args()


## Use the source dir to import helper modules
sys.path.append(args.src+'/python')
import utils
import trimadapters


## process args
TRIM_FP=args.trim_fp
TRIM_TP=args.trim_tp
SAMID = str(args.samid).split(',')
FP_ADAPTERS = str(args.fp_adapters).split(',')
TP_ADAPTERS = str(args.tp_adapters).split(',')
input = str(args.input).split(',')
output = str(args.output).split(',')


if TRIM_FP=='True' or TRIM_TP=='True':
   adapters = trimadapters.getAdapters(sam=args.sample, fp_adapters=FP_ADAPTERS, tp_adapters=TP_ADAPTERS, samid_all=SAMID)
   trimadapters.trimAdapters(infile=input, outfile=output, adapt5=adapters['fp'], adapt3=adapters['tp'], cutadapt=args.cutadapt, threads=args.threads, log=args.log)

else:
   utils.touch(output)

