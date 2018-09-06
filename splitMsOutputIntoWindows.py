#!/usr/bin/env python
import sys, gzip

msFile,numWins,winFilePrefix = sys.argv[1:]
numWins = int(numWins)

def getSnpWindowAssignments(positions,numWins):
    delta = 1.0/numWins
    if numWins > 10000:
        sys.exit("Let's not get carried away with the number of windows . . .\n")
    winStart = 0.0
    winEnd = 0+delta
    winIndex = 0
    snpWinAssignments = []
    for i in range(len(positions)):
        while not ((positions[i] > winStart or positions[i] == 0) and positions[i] <= winEnd):
            winStart = winEnd
            winEnd += delta
            winIndex += 1
            #trying to avoid some floating-point precision weirdness here
            if abs(1.0-winEnd) < 1e-9:
                winEnd = 1.0
        snpWinAssignments.append(winIndex)
    return snpWinAssignments

def getWinRange(fileIndex,numWins,delta):
    winStart = fileIndex/float(numWins)
    winEnd = winStart + delta
    if abs(1.0-winEnd) < 1e-9:
        winEnd = 1.0
    return winStart,winEnd

def getSegSitesForFiles(snpWindowAssignments,positions,numWins,outFileLs):
    segsiteCountLs = []
    segsitePositions = []
    for i in range(len(outFileLs)):
        segsiteCountLs.append(0)
        segsitePositions.append([])

    for i in range(len(snpWindowAssignments)):
        fileIndex = snpWindowAssignments[i]
        delta = 1.0/numWins
        winStart,winEnd = getWinRange(fileIndex,numWins,delta)
        windowedPosition = (positions[i]-winStart)/delta
        segsitePositions[fileIndex].append(windowedPosition)
        segsiteCountLs[fileIndex] += 1

    return segsiteCountLs,segsitePositions

def processSimulation(samples,snpWindowAssignments,positions,numWins,outFileLs):
    #first output the header information for the simulation
    segsiteCountLs,segsitePositionsLs = getSegSitesForFiles(snpWindowAssignments,positions,numWins,outFileLs)
    for i in range(len(outFileLs)):
        outFileLs[i] += "\n//\nsegsites: %s\n" %(segsiteCountLs[i])
        outFileLs[i] += "positions: " + " ".join([str(x) for x in segsitePositionsLs[i]]) + "\n"
    for sample in samples:
        for i in range(len(sample)):
            outFileLs[snpWindowAssignments[i]] += sample[i]
        for i in range(len(outFileLs)):
            outFileLs[i] += "\n"
        
if msFile == "stdin":
    isFile = False
    msStream = sys.stdin
else:
    isFile = True
    if msFile.endswith(".gz"):
        msStream = gzip.open(msFile)
    else:
        msStream = open(msFile)

header = msStream.readline()
program,numSamples,numSims = header.strip().split()[:3]
numSamples,numSims = int(numSamples),int(numSims)

#initialize list of output files
outFileLs = []
outFileNameLs = []
for i in range(numWins):
    outFileNameLs.append("%s_%s.msWin" %(winFilePrefix,i))
    outFileLs.append("./windowedMSOutput %s %s\nblah\n" %(numSamples,numSims))

processedSims = 0
#advance to first simulation
line = msStream.readline()
while not line.startswith("//"):
    line = msStream.readline()
while line:
    if not line.startswith("//"):
        sys.exit("Malformed ms-style output file: read '%s' instead of '//'. AAAARRRRGGHHH!!!!!\n" %(line.strip()))
    segsitesBlah,segsites = msStream.readline().strip().split()
    segsites = int(segsites)
    if segsitesBlah != "segsites:":
        sys.exit("Malformed ms-style output file. AAAARRRRGGHHH!!!!!\n")

    positionsLine = msStream.readline().strip().split()
    if not positionsLine[0] == "positions:":
        sys.exit("Malformed ms-style output file. AAAARRRRGGHHH!!!!!\n")
    positions = [float(x) for x in positionsLine[1:]]
    snpWindowAssignments = getSnpWindowAssignments(positions,numWins)

    samples = []
    for i in range(numSamples):
        sampleLine = msStream.readline().strip()
        if len(sampleLine) != segsites:
            sys.exit("Malformed ms-style output file %s segsites but %s columns in line: %s; line %s of %s samples AAAARRRRGGHHH!!!!!\n" %(segsites,len(sampleLine),sampleLine,i,numSamples))
        samples.append(sampleLine)
    if len(samples) != numSamples:
        raise Exception
    processSimulation(samples,snpWindowAssignments,positions,numWins,outFileLs)
    processedSims += 1
    line = msStream.readline()
    #advance to the next non-empty line or EOF
    while line and line.strip() == "":
        line = msStream.readline()
if processedSims != numSims:
    sys.exit("Malformed ms-style output file: %s of %s sims processed. AAAARRRRGGHHH!!!!!\n" %(processedSims,numSims))

for i in range(len(outFileLs)):
    outFile = open(outFileNameLs[i], "w")
    outFile.write(outFileLs[i])
    outFile.close()

if isFile:
    msStream.close()
