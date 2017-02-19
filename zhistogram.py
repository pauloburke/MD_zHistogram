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
parser = argparse.ArgumentParser(description="This software generates a histogram of densities on z axis between two groups of atoms. Created by Paulo Burke 2017")
parser.add_argument("traj", help=".xtc trajectory file.")
parser.add_argument("topology", help=".pdb topology file.")
parser.add_argument("index", help=".ndx index file.")
parser.add_argument("group1", help="Name of group 1.")
parser.add_argument("group2", help="Name of group 2.")
parser.add_argument("-b","--bins", help="Number of bins (default: numpy histogram 'auto').",type=int)
parser.add_argument("-c","--cuttoff", help="Cuttoff distance.",type=float)
parser.add_argument("-f","--frames", help="Range or number of frames to consider. Example 1-100 or 100.",type=parseRange)
parser.add_argument("-o","--output", help="File to output the histogram image and data (without extention)")
parser.add_argument("-p","--plot", help="Show plot of histogram.", action='store_true')
parser.add_argument("-s","--smooth", help="Smooth curve by applying Fourier Transform.", action='store_true', default=False)
parser.add_argument("-cv","--converge", help="Converge curve to 1.", action='store_true')
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


#setting bins
bins = "auto"
if args.bins != None:
	bins = args.bins

#get mean box size
box = np.zeros(3)
for i in xrange(len(frames)):
	box[0] += traj.unitcell_lengths[frames[i]][0]
	box[1] += traj.unitcell_lengths[frames[i]][1]
	box[2] += traj.unitcell_lengths[frames[i]][2]
	
box /= len(frames)


if args.verbose:
	print('Box mean size (nm): '+str(box))

#set output file name
if args.output != None:
	outfile = args.output
else:
	outfile = args.traj.replace(".xtc","")



#Calculation -------------------------------------------------------------------


distances = np.zeros(len(frames)*len(index[args.group1])*len(index[args.group2]))
count = 0


#calculate distances
if args.verbose:
	print('Calculating distances between '+str(args.group1)+' and '+str(args.group2)+' ...')
initfr = frames[0]
for fr in frames:
	if args.verbose:
		print(str(fr-initfr+1)+' from '+str(len(frames))+' frames',end='\r')
		sys.stdout.flush()
		
	for i in index[args.group1]:
		for j in index[args.group2]:
			distances[count] = traj.xyz[fr,i-1,2]-traj.xyz[fr,j-1,2]
			#correcting boundaries
			if distances[count]>box[2]:
				distances[count]-=box[2]
			if distances[count]<-box[2]:
				distances[count]+=box[2]
			count+=1
print()



	

if args.verbose:
	print('Calculating histogram...')
	

#setting range
if args.cuttoff != None:
	hInit = -args.cuttoff
	hEnd = args.cuttoff
else:
	hInit = np.amin(distances)
	hEnd = np.amax(distances)
	
#calculate histogram
hist,bins = np.histogram(distances,bins=bins, range=(hInit,hEnd))
#calculate bins width
dz = np.diff(bins)

#Calculating average density
rho = float(len(index[args.group1])+len(index[args.group2]))/(box[0]*box[1]*box[2])
fhist = np.zeros(len(hist))
#normalize
if args.verbose:
	print("Normalizing...")
for i in range(len(dz)):
	fhist[i] = float(hist[i])/float(sum(hist))
	volBin = dz[i]*box[0]*box[1]
	nIdeal = volBin*rho
	fhist[i] = float(fhist[i])/nIdeal
	

if args.verbose:
	print("Smoothing...")
#smooth curve with fourier
if args.smooth:
	rft = np.fft.rfft(fhist)
	rft[21:] = 0   # Note, rft.shape = 21
	fhist = np.fft.irfft(rft)

if args.verbose:
	print("Converging...")
#converge ends to 1
if args.converge:
	first = fhist[0]
	for i in range(len(fhist)):
		fhist[i]/=first	


#create plot
center = (bins[:-1]+bins[1:])/2
plt.plot(center, fhist,color='green')
plt.xlabel('nm')
plt.ylabel('g(z)')


#Save plot
if args.verbose:
	print('Saving histogram plot to '+outfile+'.png ...')
plt.savefig(outfile+".png")
 
#save data
if args.verbose:
	print('Saving histogram data to '+outfile+'.txt ...')
f = open(outfile+"txt",'w')
for i in range(len(fhist)):
	f.write(str(bins[i])+"\t"+str(fhist[i])+"\n")
f.close()

#show plot
if args.plot:
	plt.show()


if args.verbose:
	print('Done.')


