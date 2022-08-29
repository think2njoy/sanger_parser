'''
	Author: Issaac Rajan, Bioscience Core Lab @ KAUST, KSA

	Project Name: Sanger output data parser
	Start Date: Jan 2022

	Usage:
		traverse_sanger_data.py
		( Interactively provide input data path and specify output file prefix. )

	Purpose:
		To parse output data from the Sanger platform and write out the 
		data of interest to output files. Please read README.md before using this script. 

	License: 
		MIT License

	Bugs: 
		issaac.rajan@kaust.edu.sa (OR) issaac@gmail.com
'''

#!/usr/bin/env python

import sys, os, glob
from datetime import datetime

from Bio import SeqIO
from Bio.SeqIO import parse

# Procedure: (simplified)
#
# 1. Parse runlog1.txt (DONE)
# 2. From each line, identify the following: PI :: user :: sample info :: sequence length	(DONE)
# 3. Parse ab1 file. (DONE)
# 4. Print out the relevant information. (DONE)

sanger_data_path       = ''
output_filename_pos    = ''
output_filename_nonpos = ''
check_this_file  = 'check_this_file.txt' # This is for identifying files that have either or both of the keywords 'positive' or 'control' but not as 'PositiveControl'
skipped_paths    = 'skipped_paths.txt'   # This is for identifying paths that do not have some expected essential file (runlog1.txt) / folder (sequence) / etc

# make sure to remove these before you start recording what is currently missed. 
if os.path.isfile( check_this_file ): os.unlink( check_this_file )
if os.path.isfile( skipped_paths   ): os.unlink( skipped_paths   )

def parse_sanger_user_data():

   if os.path.isfile( check_this_file ): os.unlink( check_this_file )
   sanger_data_dirs = os.listdir( sanger_data_path )
#   print( 'Number of sanger data dirs: {}'.format( len( sanger_data_dirs ) ) )

   header_str = 'Date\tYearMonth\tPIID\tUserID\tRunID\tSampleName\tRL_Runlog\tRL_Fasta\tRL_Full_ab1\tRL_Trim_ab1\n'
   for output_filename in ( output_filename_pos, output_filename_nonpos):
      with open( output_filename, 'w' ) as w: w.write( header_str )

   for count, data_dir in enumerate( sanger_data_dirs, 1 ):

      if not os.path.isdir( data_dir ): pass

      full_path_for_data_dir = os.path.join( sanger_data_path, data_dir )
      full_path_for_runlog   = os.path.join( full_path_for_data_dir,  'runlog1.txt' )
      full_path_for_seq_dir  = os.path.join( full_path_for_data_dir,  'sequence' )

      if not os.path.isfile( full_path_for_runlog ):   # skip this data due to this non-existent runlog1.txt file in the respective data path. 
         mode = 'a' if os.path.isfile( skipped_paths ) else 'w'
         with open( skipped_paths, mode ) as out: out.write( full_path_for_runlog + '\n' )
         continue
      elif not os.path.isdir( full_path_for_seq_dir ): # skip this data due to this non-existent sequence/ in the respective data path. 
         mode = 'a' if os.path.isfile( skipped_paths ) else 'w'
         with open( skipped_paths, mode ) as out: out.write( full_path_for_seq_dir + '\n' )
         continue 
         
      if os.path.isfile( full_path_for_runlog ): 
         # print( full_path_for_runlog ) # For testing...
         # parse the runlog to identify: PI :: User :: Sample from each of its lines
         for line in open( full_path_for_runlog ):
            contents = line.split( '\t' )
            if not contents: continue # ignore empty lines in the runlog.
            RunID  = contents[0].strip()
            PI     = contents[2].strip()
            User   = contents[4].strip()
            Sample = contents[7].strip()
            SeqLen = contents[8].strip()
            #print( contents ) ; exit()

            try: SeqLen = int( SeqLen )
            except ValueError: continue

            if SeqLen == 0: continue

            ### BEGIN: VALID USE CASES WHERE THE PARSER NEEDS TO SKIP ###

            if ( 'positive' in Sample.lower() ) ^ ( 'control' in Sample.lower() ):
               mode = 'a' if os.path.isfile( check_this_file ) else 'w'
               mesg = 'Check this: {} --> {}'.format( RunID, Sample )
               with open( check_this_file, mode ) as out: out.write( mesg + '\n' )
               print( mesg )

            # We need something like the following:
            # date	year_month	full_path	filename	read_length	PI	User
            # Hence, parsing just the runlog1.txt will not help because we need the TIMESTAMP information. 

            #print( 'RunID :: {}, PI :: {}, User :: {}, Sample :: {}, SeqLen :: {}'.format( RunID, PI, User, Sample, SeqLen ) )

            if not Sample: continue

            if Sample.startswith( RunID + '-' + PI + '-:' ) and int( SeqLen ) == 0: continue


            Sample_ab1 = os.path.join( full_path_for_seq_dir, Sample + '.ab1' )
            Sample_seq = os.path.join( full_path_for_seq_dir, Sample + '.seq' )

            missing_files = list()

            if not os.path.isfile( Sample_ab1 ): missing_files.append( Sample_ab1 )
            if not os.path.isfile( Sample_seq ): missing_files.append( Sample_seq )

            for missing_file in missing_files:
               mode = 'a' if os.path.isfile( skipped_paths ) else 'w'
               with open( skipped_paths, mode ) as out: out.write( missing_file + '\n' )
               
            if missing_files: continue

            ### END: VALID USE CASES THAT REQUIRE A SKIP ###

            # compute read length from Sample_seq file. Cannot use biopython module bcos these single-seq fasta files are without a sequence header. That's why.
            SeqLen_from_Fasta = 0
            with open( Sample_seq ) as fasta_fp:
               for line in fasta_fp:
                  if line.startswith( '>' ): continue # takes care of situations where the single sequence file is devoid of a header
                  SeqLen_from_Fasta += len( line[:-1] )

            rmode_RLs = { 'abi':0, 'abi-trim':0 }
            for rmode in ( 'abi', 'abi-trim' ):
               with open(Sample_ab1, 'rb') as handle:
                  record = next( SeqIO.parse( handle, rmode ) )
                  sequence = record.__dict__['_seq']._data
                  read_length = len(sequence)
                  rmode_RLs[ rmode ] = read_length
                  time_mod_time_since_epoch = os.path.getmtime( Sample_ab1 )
                  datetime_time = datetime.fromtimestamp( time_mod_time_since_epoch )
                  user_friendly_fname = os.path.basename( Sample_ab1 ).split('.')[0]

            output_filename = output_filename_pos if 'PositiveControl' in Sample else output_filename_nonpos
            wmode = 'a' if os.path.isfile( output_filename ) else 'w'
            with open( output_filename, wmode ) as w: w.write( '{0}-{1}-{2}\t{0}-{1}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n'.format( datetime_time.year, str(datetime_time.month).zfill(2), str(datetime_time.day).zfill(2), PI, User, RunID, user_friendly_fname, SeqLen, SeqLen_from_Fasta, rmode_RLs['abi'], rmode_RLs['abi-trim'] ) )
    
if __name__ == '__main__':

   sanger_data_path = input( 'Specify sanger data path (needed): ' )
   output_filename_prefix  = input( 'Specify the prefix for sanger output file names (optional): ' )

   output_filename_pos    = output_filename_prefix + '_pos.txt'
   output_filename_nonpos = output_filename_prefix + '_nonpos.txt'

   for output_filename in ( output_filename_pos, output_filename_nonpos ):
      if os.path.isfile( output_filename ):
         print( 'A file by the name {} found. Do you want to overwrite this?'.format( output_filename ) )
         response      = input( 'Okay to overwrite [ yes / y / no / n ] --> ' )
         if response.lower() not in ( 'y', 'yes' ): sys.exit( 'Quitting the program based on your response {}'.format( response ) )

   parse_sanger_user_data()
