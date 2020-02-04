#!/usr/bin/python
##Takes the output of gromacs contact analysis and puts it into a "friendly" format.
##Usage: ContactFreq.py gromacsOutput.xpm systemSize [framestride]
##Inputs:
##	gromacsOutput.xpm: output from using -frames with mdmat 
##	systemSize: number of atoms in the contact matrix
##	framestride: (optional) only read in every n frame
##Outputs:
##	FrequencyMap.dat: text file containing the frequency of each contact. Contains systemSize lines of length systemSize; contact 1->1 is line 1, number 1, contact 2->1 is line 2, number 1, etc.
##	Matplotlib image: representation of the frequency matrix shown through matplotlib.pyplot.imshow(). Number below tol (hardcoded) are displayed as white.
##Written by Joseph Clayton (jclayto1@hawk.iit.edu), 2019
##Dec, 2019: Fixed a bug where the first line was ignored when reading, leading to a one-residue shift in the matrix

import sys
import numpy as np
import matplotlib.pyplot as plt

##Hardcoded parameters
cutoff = 0.4		#nm; does not have to be the same used in mdmat command
tol = 1e-2		#Hide frequencies below this value

#Check and prepare inputs
if(len(sys.argv) < 3):
	print('Usage: %s gromacsOutput.xpm systemSize [frameStride]'%sys.argv[0])
	sys.exit() 

try:
	numAtoms = int(sys.argv[2])
except:
	print('Error: unable to parse systemSize parameter')
	sys.exit()

striding = False
try:
	frameStride = int(sys.argv[3]) -1
	striding = True
except:
	pass

#Initialize matrix
contactMap = np.zeros((numAtoms,numAtoms))

#Open the file
try:
	contactFile = open(sys.argv[1],'r')
except:
	print('Error: unable to open %s'%sys.argv[1])
	sys.exit()

#Read in the scale
contactChars = []
headingSize = 0
for line in contactFile:
	#Skip heading information (GROMACS command, etc):
	headingSize+=1
	if('x-axis' in line):	break		#Ends when the contact scale has been read fully
	if (line[0] == '/'): continue		#Comments
	if (line.split()[0] == 'static'): continue #Definition of scale
	if (line.split()[1].isdigit()): continue#Definition of matrix
	#If a color corresponds to a close distance, save it 
	distance=float(line.split()[5].replace('"',''))
	if distance<cutoff:
		contactChars.append(line[1])
#Continue reading the heading
totalFrame= 0
j=0
skiped=0
for line in contactFile:
	if ('axis' in line):
		headingSize+=1
	else:
		#Loop through minimum distances
		for char in ['"',',','\n']: line=line.replace(char,'')
		if (len(line)!=numAtoms): continue
		for i,dist in enumerate(line):
			if (dist in contactChars): contactMap[i,j]+=1
		j+=1
		if (j==numAtoms): 
			totalFrame+=1
			j=0
			if (striding): 
				for i in xrange((numAtoms+headingSize)*frameStride): 
					try:
						next(contactFile)
					except:
						break	#EOF
#Calculate frequency from contact map
freqMap = contactMap[:,::-1] / totalFrame
np.savetxt('FrequencyMap.dat',freqMap,delimiter=' ',fmt='%1.4f')

#Visualize the matrix
freqMap[freqMap<tol] = np.nan
plt.imshow(freqMap,origin='lower',aspect='auto')
plt.colorbar()
plt.show()
