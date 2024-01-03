#############################################################################################################
# pipeline_template
#
# This program is free software that contains third party software subject to various licenses, 
# namely, the GNU General Public License version 3 (or later), the GNU Affero General Public License 
# version 3 (or later), and the LaTeX Project Public License v.1.3(c). A list of the software contained 
# in this program, including the applicable licenses, can be accessed here: 
#
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# Program:  run-report.r
# Version:  0.1
# Author:   Sami R. Cherikh
# Purpose:  Renders summary report for pipeline
# Input:    <COMMAND LINE INPUTS>
# Output:
#############################################################################################################

args = commandArgs(trailingOnly = TRUE)
if(length(args)>0) {
    rmd.file = args[1]
    predir = args[2]
    outdir = args[3]
}else{
    stop("ERROR: No rmd file, predir, and/or outdir supplied.")
}
# cp worlkflow png from srcdir to predir
srcdir = paste0(dirname(rmd.file),'/..')
system(paste0("cp ",srcdir,"/../doc/*_workflow.png ",predir))
# knit rmarkdown html report
rmarkdown::render(rmd.file, output_format="slidy_presentation", output_dir=outdir, params = list(predir = predir))
