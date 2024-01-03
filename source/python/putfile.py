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
# Program:  putfile.py
# Version:  0.1
# Author:   Sami R. Cherikh, Travis L. Jensen
# Purpose:  Interface for uploading files to the cloud
# Input:    N/A
# Output:   N/A
#############################################################################################################

import subprocess
import logging
import utils


## Upload to an Amazon S3 or Google Cloud bucket, preserving folder structure
def cloudUpload(file, destination, prog, doencrypt='False', rm=True):
    if destination[-1] != '/':
        destination += '/'

    if not isinstance(file, list):
        file = file.split(' ')

    for f in file:
        if('aws' in prog):
            cmd = '%s s3 cp %s %s' % (prog, f, destination+f)
        else:
            cmd = '%s cp %s %s' % (prog, f, destination+f)
        ## command to remove encrypted file
        if (doencrypt=='True' and rm):
            cmd_rm = 'rm %s' % (f)

        ## Attempt to upload
        try:
            ## upload
            utils.logging_call(cmd, shell=True)
            ## remove encrypted file after uploading
            if (doencrypt=='True' and rm):
                utils.logging_call(cmd_rm, shell=True)
        except subprocess.CalledProcessError:
            logging.error('Cloud upload failed. See above for more details.')
            exit(1)
    return(True)



## Use different functions depending on cloud provider (AWS and Google)
def upload(file, destination, prog, doencrypt='False', rm=True):
    #upload_func = {cloud: cloudUpload}
    #upload_func[cloud](file=file, destination=destination, prog=prog, rm=rm)
    cloudUpload(file=file, destination=destination, prog=prog, doencrypt=doencrypt, rm=rm)