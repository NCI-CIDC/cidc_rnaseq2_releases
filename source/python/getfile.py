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
# Program:  getfile.py
# Version:  0.1
# Author:   Sami R. Cherikh, Travis L. Jensen
# Purpose:  Interface for pulling input files
# Input:    N/A
# Output:   N/A
#############################################################################################################

import os
import subprocess
import shutil

import utils


## Download file(s) from AWS S3 or Google Cloud storage, return path to result
def downloadCloudFile(in_file,out_file,cloud):
	res = []
	for f in in_file:
		if('aws' in cloud):
			#cmd = cloud + ' s3 cp ' + f + ' ' + os.getcwd() + '/' + os.path.basename(f)
                        cmd = cloud + ' s3 cp ' + f + ' ' + out_file
		else:
			#cmd = cloud + ' cp ' + f + ' ' + os.getcwd() + '/' + os.path.basename(f)
                        cmd = cloud + ' cp ' + f + ' ' + out_file
		utils.logging_call(cmd, shell=True)
		res.append(os.path.basename(f))
	return(res)



## Construct SRA URL (see https://www.ncbi.nlm.nih.gov/books/NBK158899/)
def getSRAUrl(accession):
	res = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/'+accession[0:6]+'/'+accession+'/'+accession+'.sra'
	return(res)

	
## Download file(s) from the sequence read archive, return path to result
def downloadSRAFile(accession, fastq_dump, end=1):
	res = []
	for f in accession:
	
		## Download from the NCBI ftp site, extract fastq.gz, remove downloaded file
		## pefetch is more reliable than fastq-dump to retrieve the sra file
		if(end==1):
			cmd = 'prefetch -v ' + f + ' -O ' + os.getcwd() + ' && ' + fastq_dump + ' --sra-id ' + os.getcwd() + '/' + f + '/' + f + '.sra --gzip --outdir ' + os.getcwd() + ' --threads 8 && rm -r ' + os.getcwd() + '/' + f
		else:
			cmd = 'prefetch -v ' + f + ' -O ' + os.getcwd() + ' && ' + fastq_dump + ' --sra-id ' + os.getcwd() + '/' + f + '/' + f + '.sra --gzip --outdir ' + os.getcwd() + ' --threads 8 --split-3 && rm -r ' + os.getcwd() + '/' + f

		## If download times out on the first try, give it up to 2 more tries
		for i in range(3):
			try:
				utils.logging_call(cmd, shell=True)
			except subprocess.CalledProcessError as e:
				print('Timeout: ' + e.cmd)
				if (i == 2):
					raise e
			else:
				break
		## Handle paired end output files from splitting accession if paired end
		if(end==1):
			res.append(os.getcwd() + '/' +f+'.fastq.gz')
		else:
			res.append(os.getcwd() + '/' + f + '_1.fastq.gz')
			res.append(os.getcwd() + '/' + f + '_2.fastq.gz')
	return(res)


def getFile(in_file, out_file, cloud_prog, fastq_dump_prog, openssl_prog, pw, hash, end=1):
	## Set flags
	file_downloaded = False
	file_decrypted  = False
	file_sra = False
	
	## Convert to list so we can handle multiple files
	in_file = in_file.split(';')
	
	## Download file (if needed)
	if (in_file[0][0:5] == 's3://' or in_file[0][0:5] == 'gs://'):
		in_file = downloadCloudFile(in_file, out_file, cloud_prog)
		file_downloaded = True
	elif(in_file[0][0:3] == 'SRR'):
		in_file = downloadSRAFile(in_file, fastq_dump_prog, end)
		file_downloaded = True
		file_sra = True


	## Decrypt file (if needed)
	if (in_file[0][-4:] == '.enc'):
		in_file = utils.decryptFile(in_file, openssl_prog, pw, hash)
		file_decrypted = True
	
	## If we downloaded & decrypted, remove the encrypted file(s)
	if (file_downloaded and file_decrypted):
		[os.remove(s+'.enc') for s in in_file]

	# ## Merge if multiple downloaded files, only keep merged file
	# if (len(in_file) > 1 and file_downloaded):
	# 	utils.cat(in_file, out_file)
	# 	[os.remove(s) for s in in_file]

	## Merge if multiple downloaded files, only keep merged file
	if (len(in_file) > 1 and file_downloaded and not file_sra):
		utils.cat(in_file, out_file)
		[os.remove(s) for s in in_file]

	## Handle copying split PE SRA files
	elif (len(in_file) > 1 and file_downloaded and file_sra):
		shutil.copy2(in_file[0], out_file)
		shutil.copy2(in_file[1], out_file.replace('_1.fastq.gz', '_2.fastq.gz'))
		[os.remove(s) for s in in_file]
	
	## Merge if multiple local files, keep unmerged files
	elif (len(in_file) > 1 and not file_downloaded):
		utils.cat(in_file, out_file)
	
	## Move if single downloaded file
	elif (len(in_file) == 1 and file_downloaded):
		#os.rename(in_file[0], out_file)
                print("Download successful")	
	## Make a copy if single local file
	elif (len(in_file) == 1 and not file_downloaded):
		#cmd = 'cp '+in_file[0]+' '+out_file
		#subprocess.run(cmd, shell=True, check=True)
		shutil.copy2(in_file[0], out_file)
		
	
	## If everything runs successfully, return 0 
	return(0)
