3D-DART server README file.

Thank you for choosing the 3D-DART server for modeling your DNA structures. 
You have successfully extracted the archive with your results and are now
faced with a set of directories and a few files. 
The modeling job that was performed for you consists out of different
scripts linked together and executed in batch-mode. Every script in the 
batch saves its output in its own directory that bears the number of the 
job and the name of the script. A brief description of the directories:

jobnr1-FileSelector:    	Contains the files you uploaded if you did so,
					 		otherwise it remains empty. 
jobnr2-BuildNucleicAcids: 	Generates a canonical starting structure of
						  	the DNA with the desired sequence
jobnr3-X3DNAAnalyze: 		Runs a 3DNA analysis over the structure generated in
					 		job 2 to obtain initial starting parameters for modeling
jobnr4-ModelNucleicAcids: 	The actual modeling step. Models are represented as
						  	base-pair(step) parameter files.
jobnr5-BuildNucleicAcids: 	Converts the parameter files from job 4 into PDB
 						  	structure files.
jobnr6-X3DNAanalyze: 		Runs a 3DNA analysis over the modeled structures. 
jobnr7-NABendAnalyze: 		Runs a global bend analysis over the modeled structures.
jobnr8-PDBeditor: 			Makes some last-minute changes to the PDB files. 
							This directory contains the final structures.

In the main directory you will find the log file of your job named 'dart.out'
In case you encounter unexpected errors please look into the log file.
If the source of the problem remains unclear pleas e-mail the log file to
the maintainer of the 3D-DART web-server.

If you use 3D-DART please cite:

	M. van Dijk and A.M.J.J. Bonvin (2009) "3D-DART: a DNA structure modelling server", Nucl. Acids Res.
    37 (Web Server Issue):W235-W239, doi:10.1093/nar/gkp287

Kind Regards,

The 3D-DART server team
  