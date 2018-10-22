# parallel_rosetta

This code takes a pre-installed copy of Rosetta Ab Initio and runs it in parallel.
You must provide paths to the fragment files, fasta, and PsiPred2 files, and
you must have working copies of "AbInitioRelax" and "extract_pdbs". If you clone
my "IDP_analysis" class then you can easily analyze all of your structures
using this script as well.

Requirements: 	Rosetta
		Python 2.7
		Joblib
		Numpy

I find Joblib to be a little difficult to install sometimes, so I think it's easiest
to download a copy of Miniconda 3.x, and then activate the environment with *activate py2*
This will give you a python 2.7 sandbox and it will be easy to install all of the 
requirements with pip.
