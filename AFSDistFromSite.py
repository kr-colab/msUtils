#!/usr/bin/env python
import sys

msFile, targetPos, maxDistanceFromTarget = sys.argv[1:]
targetPos = float(targetPos)
maxDistanceFromTarget = float(maxDistanceFromTarget)

def processSimulation(samples, positions, targetPos, maxDistanceFromTarget, afs):
    freqH = {}
    for i in range(len(samples[0])):
        if abs(targetPos - positions[i]) < maxDistanceFromTarget:
            freqH[i] = {}
            for sample in samples:
                if not freqH[i].has_key(sample[i]):
                    freqH[i][sample[i]] = 0
                freqH[i][sample[i]] += 1
    for i in freqH.keys():
        alleles = sorted(freqH[i].keys())
        if alleles == ['0', '1']:
            freq = freqH[i]['1']
            afs[freq] += 1

if msFile == "stdin":
    isFile = False
    msStream = sys.stdin
else:
    isFile = True
    msStream = open(msFile)

header = msStream.readline()
program,numSamples,numSims = header.strip().split()[:3]
numSamples, numSims = int(numSamples), int(numSims)

processedSims = 0
#advance to first simulation
line = msStream.readline()
while line.strip() != "//":
    line = msStream.readline()
afs = -1
while line:
    if line.strip() != "//":
        sys.exit("Malformed ms-style output file: read '%s' instead of '//'. AAAARRRRGGHHH!!!!!\n" %(line.strip()))
    segsitesBlah,segsites = msStream.readline().strip().split()
    segsites = int(segsites)
    if segsitesBlah != "segsites:":
        sys.exit("Malformed ms-style output file. AAAARRRRGGHHH!!!!!\n")

    if segsites != 0:
        positionsLine = msStream.readline().strip().split()
        if not positionsLine[0] == "positions:":
            sys.exit("Malformed ms-style output file. AAAARRRRGGHHH!!!!!\n")
        positions = [float(x) for x in positionsLine[1:]]

        samples = []
        for i in range(numSamples):
            sampleLine = msStream.readline().strip()
            if len(sampleLine) != segsites:
                sys.exit("Malformed ms-style output file %s segsites but %s columns in line: %s; line %s of %s samples AAAARRRRGGHHH!!!!!\n" %(segsites, len(sampleLine), sampleLine, i, numSamples))
            samples.append(sampleLine)
        if len(samples) != numSamples:
            raise Exception
        if afs == -1:
            afs = [0]*len(samples)
        processSimulation(samples, positions, targetPos, maxDistanceFromTarget, afs)
    processedSims += 1
    line = msStream.readline()
    #advance to the next non-empty line or EOF
    while line and line.strip() == "":
        line = msStream.readline()
#print afs[1:]
denom = float(sum(afs[1:]))
for i in range(1,len(afs)):
    print "%i %le" %(i, afs[i]/denom)
if processedSims != numSims:
    sys.exit("Malformed ms-style output file: %s of %s sims processed. AAAARRRRGGHHH!!!!!\n" %(processedSims, numSims))

if isFile:
    msStream.close()
