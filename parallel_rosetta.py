import os
import subprocess
import numpy as np
import shutil
import random
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory

global total_structures ; total_structures = 50000
global nruns            ; nruns = 10
global nstrucs          ; nstrucs = int(float(total_structures)/nruns)

mode = ['RUN'] # RUN, EXTRACT, ANALYZE

def run_abinitio(d):
	if d < 10:
		dir = os.getcwd()+"/0"+str(d)+'/'
	else:
		dir = os.getcwd()+"/"+str(d)+'/'
	try:
		os.mkdir(dir)
	except:
		print dir, "already exists!"
	os.chdir(dir)
	shutil.copy2('/home/dillion/Dropbox/Reflectin/structure_prediction/rosetta/2_8/00/rob_motif2_8.fasta',               dir+'rob_motif2_8.fasta')
	shutil.copy2('/home/dillion/Dropbox/Reflectin/structure_prediction/rosetta/2_8/00/rob_motif2_8.frag-03_05.200_v1_3', dir+'rob_motif2_8.frag-03_05.200_v1_3')
	shutil.copy2('/home/dillion/Dropbox/Reflectin/structure_prediction/rosetta/2_8/00/rob_motif2_8.frag-09_05.200_v1_3', dir+'rob_motif2_8.frag-09_05.200_v1_3')
	shutil.copy2('/home/dillion/Dropbox/Reflectin/structure_prediction/rosetta/2_8/00/rob_motif2_8.psipred_ss2',         dir+'rob_motif2_8.psipred_ss2')
	rand = random.randint(1,10000000)

	FNULL = open(os.devnull, 'w')
	subprocess.call(['AbinitioRelax', \
		'-in::file::fasta rob_motif2_8.fasta', \
		'-in:file:frag3 rob_motif2_8.frag-03_05.200_v1_3', \
		'-in:file:frag9 rob_motif2_8.frag-09_05.200_v1_3', \
		'-abinitio:relax', \
		'-constant_seed', \
		'-jran ', str(rand), \
		'-relax:fast', \
		'-use_filters true', \
		'-psipred_ss2 rob_motif2_8.psipred_ss2', \
		'-nstruct ', str(nstrucs), \
		'-out:file:silent fast_relax_1.out'], stdout=FNULL, stderr=subprocess.STDOUT)

	return 0

def run_extract(n):
	dir = os.getcwd()+"/0"+str(n)+'/'
	try:
		os.mkdir(dir+'pdbs')
	except:
		print "pdbs/ already exists"
	FNULL = open(os.devnull, 'w')
	subprocess.call(['extract_pdbs', '-in::file::silent', dir+'fast_relax_1.out', '-out::prefix', dir+'pdbs/ref'], stdout=FNULL, stderr=subprocess.STDOUT)
	for i in range(1,5001):
		j = i+n*5000
		if i < 10:
			n1 = '000'+str(i)+'.pdb'
			n2 = '000'+str(j)+'.pdb'
		elif i < 100:
			n1 = '00'+str(i)+'.pdb'
			n2 = '00'+str(j)+'.pdb'
		elif i < 1000:
			n1 = '0'+str(i)+'.pdb'
			n2 = '0'+str(j)+'.pdb'
		elif i < 10000:
			n1 = str(i)+'.pdb'
			n2 = str(j)+'.pdb'

		try:
			shutil.copy2(dir+'pdbs/refS_0000'+n1,'/home/dillion/Dropbox/Reflectin/structure_prediction/rosetta/2_8/pdbs/refS_0000'+n2)
		except:
			print dir+'pdbs/refS_0000'+n1, "does not exist"

def combine_scorefiles():
	pass

def analyze(n):
	pass 

if __name__ == "__main__":

	if 'RUN' in mode:
		runs = Parallel(n_jobs=nruns)(delayed(run_abinitio,has_shareable_memory)(d) for d in range(0,nruns))

	elif 'EXTRACT' in mode:
		runs = Parallel(n_jobs=nruns)(delayed(run_extract,has_shareable_memory)(d) for d in range(0,nruns))
		combine_scorefiles()

	elif 'ANALYZE' in mode:
		runs = Parallel(n_jobs=nruns)(delayed(analyze,has_shareable_memory)(d) for d in range(0,nruns))

	else:
		print "keyward", mode, "not recognized"
