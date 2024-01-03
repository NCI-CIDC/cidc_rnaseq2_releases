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
# Program:  init-encrypt-archive.py
# Version:  0.1
# Author:   Sami R. Cherikh, Travis L. Jensen
# Purpose:  Initializer script for encrpyting and archiving output
# Input:    N/A
# Output:   N/A
#############################################################################################################

## Import modules
import sys, argparse

parser=argparse.ArgumentParser()

parser.add_argument("--src", help="Source dir")
parser.add_argument("--doencrypt", help="Do encrypt or not")
parser.add_argument("--openssl", help="Encryption program", default='', required=False, nargs='?', const='')
parser.add_argument("--encrypt-pass", help="Encryption pswd", default='', required=False, nargs='?', const='')
parser.add_argument("--hash", help="Hash alg", default='', required=False, nargs='?', const='')
parser.add_argument("--doarchive", help="Do archive or not")
parser.add_argument("--archive", help="Cloud destination", default='', required=False, nargs='?', const='')
parser.add_argument("--cloud", help="Cloud program")
parser.add_argument("--output", help="Output")

args=parser.parse_args()

## Use the source dir to import helper modules
sys.path.append(args.src+'/python')
import utils
import putfile


output = str(args.output).split(',')


## Encrypt and/or upload if necessary
if args.doencrypt=='True': output = utils.encryptFile(file=utils.toList(output), openssl=args.openssl, password=args.encrypt_pass, hash=args.hash)
if args.doarchive=='True': putfile.upload(file=utils.toList(output), destination=args.archive, prog=args.cloud, doencrypt=args.doencrypt)