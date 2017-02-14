#!/usr/bin/python

from __future__ import print_function
import sys,argparse
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import re

def read_index(index_file):
	assert (index_file.split('.')[1]=='ndx'),"Wrong file format. Index file shoud have a .ndx extention."
	index = {}
	f = open(index_file,'r')
	key = ''
	for line in f.readlines():
		if line[0] == '[':
			key = line[line.find("[")+2:line.find("]")-1]
			index[key] = []
		else:
			for j in [int(i) for i in line.split(' ') if i<>'' and i<>'\n']:
				index[key].append(j)
	return index
	
def parseRange(string):
    m = re.match(r'(\d+)(?:-(\d+))?$', string)
    # ^ (or use .split('-'). anyway you like.)
    if not m:
        raise argparse.ArgumentTypeError("'" + string + "' is not a range or number. Expected forms like '0-5' or '2'.")
    start = m.group(1)
    end = m.group(2) or start
    return list(range(int(start,10), int(end,10)+1))


#Reading and treating arguments ----------------------------------

import argparse
parser = argparse.ArgumentParser(description="This software generates a histogram of densities on z axis between two groups of atoms.")
parser.add_argument("traj", help=".xtc trajectory file.")
parser.add_argument("topology", help=".pdb topology file.")
parser.add_argument("index", help=".ndx index file.")
parser.add_argument("group1", help="Name of group 1.")
parser.add_argument("group2", help="Name of group 2.")
parser.add_argument("-f","--frames", help="Range or number of frames to consider. Example 1-100 or 100.",type=parseRange)
parser.add_argument("-c","--cuttoff", help="Cuttoff distance.",type=float)
parser.add_argument("-s","--show", help="Show plot of histogram.", action='store_true')
parser.add_argument("-o","--output", help="File to output the histogram (.png)")
parser.add_argument("-v","--verbose", help="Make the software more verbose.",action="store_true")
args = parser.parse_args()



# Treating arguments ----------------------------------------------------------------

#reading trajectories
if args.verbose:
	print('Reading trajectory from '+args.traj+' and topology from '+args.topology+'...')
traj = md.load(args.traj,top=args.topology)
if args.verbose:
	print('Trajecroty data:')
	print('# frames: '+str(traj.n_frames))
	print('# atoms: '+str(traj.n_atoms))
	print('# residues: '+str(traj.n_residues))
	

#setting frames for calculation
frames = []
if args.frames == None:
	frames = range(traj.n_frames)
	if args.verbose:
		print("Considering all frames for calculation.")
else:
	
	if len(args.frames) == 1:
		assert (args.frames[0]<=traj.n_frames),"The given max frame is bigger than number of frames in file ("+str(traj.n_frames)+")"
		frames = range(args.frames[0])
		if args.verbose:
			print("Considering frames from 0 to "+str(args.frames[0]-1)+" for calculation.")
	else:
		assert (max(args.frames)<traj.n_frames),"The given range of frames is bigger than number of frames in file ("+str(traj.n_frames)+")"
		frames = args.frames


#reading indexes
if args.verbose:
	print('Reading atoms indexes from '+args.index+'...')
index = read_index(args.index)

assert (args.group1 in index.keys()),"Group "+args.group1+" not found in index file. Options are: "+str(index.keys())
assert (args.group2 in index.keys()),"Group "+args.group2+" not found in index file. Options are: "+str(index.keys())

#setting cuttof
cut = 0
if args.cuttoff != None:
	cut = args.cuttoff
print(cut)

#get box size
box = traj.unitcell_lengths[0]
if args.verbose:
	print('Box size (nm): '+str(box))

#set output file name
if args.output != None:
	outfile = args.output
else:
	outfile = args.traj.replace(".xtc",".png")



#Calculation -------------------------------------------------------------------


distances = np.zeros(len(frames)*len(index[args.group1])*len(index[args.group2]))
count = 0

if args.verbose:
	print('Calculating distances...')
for fr in frames:
	if args.verbose:
		print(str(fr+1)+' from '+str(len(frames))+' frames',end='\r')
		sys.stdout.flush()
		
	for i in index[args.group1]:
		for j in index[args.group2]:
			distances[count] = traj.xyz[fr,i-1,2]-traj.xyz[fr,j-1,2]
			count+=1
print()


if cut:
	if args.verbose:
		print('Applying cuttoff distance...')
	distances = distances[distances<=cut]
	distances = distances[distances>=-cut]
	

if args.verbose:
	print('Calculating histogram...')
	

hist,bins = np.histogram(distances, bins='auto')
dz = np.diff(bins)
#Calculating average density
rho = float(traj.n_atoms)/(box[0]*box[1]*box[2])

#normalize
for i in range(len(dz)):
	hist[i] = hist[i]/len(frames)
	volBin = dz[i]*box[0]*box[1]
	nIdeal = volBin*rho
	hist[i] = hist[i]/nIdeal


center = (bins[:-1]+bins[1:])/2
width = 0.7*(bins[1]-bins[0])
plt.bar(center, hist, align = 'center', width = width,facecolor='green', edgecolor='green', alpha=0.5)
plt.xlabel('nm')
plt.ylabel('N')



if args.verbose:
	print('Saving histogram plot to '+outfile+'...')
plt.savefig(outfile) 


if args.show:
	plt.show()


if args.verbose:
	print('Done.')


