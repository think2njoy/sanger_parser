
# Sanger Data Output Parser


## Introduction

The purpose of this script is to be able to parse sanger data in order to obtain key parameters such as  
the time of sequencing, information about user groups, sequencing outcomes in terms of read length, etc.  
Such information when plotted graphical can lend crucial insights about the overall sequencing performance  
over time, facilitate comparsion of data quality among different user groups, etc. Further, this script  
neatly separates positive control and non-positive control output data.  


## [I] INPUT DATA

This script assumes the following structure for the SANGER date folder:  
  
```
<PATH_TO>/<DIR:SANGER_input_folder>/

   <DIR:RunID1>/
      "runlog1.txt"
      <DIR:"sequence">/
         sample1.ab1
         sample1.seq
         ...
   ... ( There can be MANY such RunIDx folders under <DIR:SANGER_input_folder> )
   <DIR:RunIDn>/ folders...
```
  
**Note:**  

1. Names under quotes imply that such a entity actually named so is expected to exist in the required data path structure.  

For example, there should be a runlog1.txt file under each directory corresponding to data from one of the runs. 
Similarly, there should be a subdirectory named sequence/ under each directory corresponding to data from one of the runs. 

2. Other unquoted names like SANGER_input_folder, RunID1, RunIDn are just place holders and are not actual expected folder names. 

3. Parsing of runlog1.txt file is done by this script, for the following information: 

RunID from the 1st column, PIID from the 3rd column, UserID from the 5th column, Sample name from the 8th column and  
sample sequence length from the 9th column. 

4. Basically, runlog1.txt serves to maintain some metadata for each RunID, which can be used to lookup for information   
as done in this case. Current version of this script requires this information to exist. 

## [2] OUTPUT DATA

Output (if any) is written to file(s). Files are named based on the prefix specified. Positive control data files have the 
suffix _pos.txt, whereas non-positive control data output file has the suffix _nonpos.txt. These output files have the 
following data columns. 

Date	YearMonth	PIID	UserID	RunID	SampleName	RL_Runlog	RL_Fasta	RL_Full_ab1	RL_Trim_ab1 

**Note:**

RL_Runlog stands for read length obtained by parsing the file runlog1.txt. 
RL_Fasta stands for read length obtained by parsing the sample fasta (.seq) sequence file. 
RL_Full_ab1 stands for read length obtained by parsing the sample ab1 (.ab1) file. 
RL_Trim_ab1 stands for trimmed read length obtained by parsing the sample ab1 (.ab1) file.

Other outputs generated are as follows:

_check_this_file.txt_ - This file lists skipped items due to keyword constraints specified in the program. 
_skipped_paths.txt_   - This file lists skipped files due to incomplete / missing data. 


## [3] SOFTWARE REQUIREMENTS

For this script to run it requires biopython to be available. Specifically, script was tested under biopython v1.79 and python v3.10.4. 


## [4] LICENSE  

This software is released under MIT license.  
  

## [5] CONTACT  
  
For reporting bugs contact: issaac.rajan@kaust.edu.sa (OR) issaac@gmail.com  

