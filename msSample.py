import sys, gzip, random
import numpy as np

class msSample:
	def __init__(self, segsites, pos,haplos):
		self.segsites = segsites
		self.positions = pos
		self.haplotypes = haplos
		self.n = len(self.haplotypes)
	
	def convertHaplosFiniteSites(self,seqlength):
		newH = np.zeros((self.n,seqlength))
		if self.segsites > 0:
			newPos = [np.floor(x*seqlength) for x in self.positions]
			for i in list(range(0,self.n)):
				for j in list(range(0,self.segsites)):
					newH[i,int(newPos[j])]=self.haplotypes[i][j]
				
		self.haplotypes=newH		
			

class msFile:
	def __init__(self):
		self.repNumber=0
		self.samples = []
		
	def readFile(aFileName):
		sampleList=[]
		msStream = open(aFileName)
		header = msStream.readline().strip().split()
		program,numSamples,numSims = header[:3]
		if len(header) > 3:
			otherParams = " " + " ".join(header[3:])
		else:
			otherParams = ""
		
		numSamples, numSims = int(numSamples),int(numSims)
		#advance to first simulation
		line = msStream.readline()
		while line.strip() != "//":
			line = msStream.readline()
		repLs = []
		while line:
			if line.strip() != "//":
				sys.exit("Malformed ms-style output file: read '%s' instead of '//'. AAAARRRRGGHHH!!!!!\n" %(line.strip()))
			s = int(msStream.readline().strip().split()[-1]) #segsites line
			h=[]
			positionsLine=[]
			if s > 0:
				positionsLine = msStream.readline().strip().split()
				if not positionsLine[0] == ("positions:"):
					sys.exit("Malformed ms-style output file. AAAARRRRGGHHH!!!!!\n")
					positionsLine.pop(0)
					positionsLine = [float(x) for x in positionsLine]
		
			
				for i in range(numSamples):
					currLine = msStream.readline()
					h.append([int(x) for x in list(currLine.strip())])
			line = msStream.readline()
			#advance to the next non-empty line or EOF
			while line and line.strip() == "":
				line = msStream.readline()
			sampleList.append(msSample(s,positionsLine,h))	
			
		msStream.close()
		return sampleList
			
