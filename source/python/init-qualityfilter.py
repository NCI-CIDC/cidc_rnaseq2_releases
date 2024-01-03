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
# Program:  init-qualityfilter.py
# Version:  0.1
# Author:   Sami R. Cherikh, Travis L. Jensen
# Purpose:  Initializer script quality trimming input reads based on quality cutoff threshold
# Input:    N/A
# Output:   N/A
#############################################################################################################

## Import modules
import os
import sys, argparse

parser=argparse.ArgumentParser()

parser.add_argument("--src", help="Source dir")
parser.add_argument("--ends", help="Single or paired end indication")
parser.add_argument("--threads", help="Threads")
parser.add_argument("--trimmomatic", help="Quality filter prog")
parser.add_argument("--qual-cutoff", help="Base quality filter cutoff")
parser.add_argument("--input", help="Input")
parser.add_argument("--output", help="Output")
parser.add_argument("--log", help="Log file")

args=parser.parse_args()


## Use the source dir to import helper modules
sys.path.append(args.src+'/python')
import utils


## process args
ENDS = str(args.ends).split(',')
threads = args.threads
QUAL_CUTOFF = args.qual_cutoff
input = str(args.input).split(',')
output = str(args.output).split(',')


# Build command (if quality cutoff is 0 just copy files to their new locations)
if (int(QUAL_CUTOFF) > 0):
   paired_end_flag = 'PE' if len(ENDS) == 2 else 'SE'
   cmd_input = ' '.join(utils.toList(input))
   cmd_output = ' '.join(utils.toList(output))
   cmd = '%s %s -threads %s %s %s SLIDINGWINDOW:4:%s' % (args.trimmomatic, paired_end_flag, threads, cmd_input, cmd_output, QUAL_CUTOFF)
else:
   ## Join together all copy commands
   cmd = ['cp %s %s' % (x,y) for x,y in zip(utils.toList(input), utils.toList(output))]
   cmd = '; '.join(cmd)

## Run command
print(cmd)
os.system("echo '"+cmd+"' >> "+args.log)
utils.logging_call(cmd, shell=True)
