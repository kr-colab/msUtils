#!/usr/bin/env python
import sys, gzip, random

msFileName1, numReps1, msFileName2, numReps2, shuffle = sys.argv[1:]
numReps1, numReps2 = int(numReps1), int(numReps2)
if not shuffle in ["shuffle", "no_shuffle"]:
    sys.exit("shuffle must be set to either 'shuffle' or 'no_shuffle'. AAAARRRRGGGGHHHHHHHHHH!!!!\n")

def readAllMSRepsFromFile(msFileName):
    msStream = open(msFileName)

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
        repStr = ["\n//"]
        repStr.append(msStream.readline().strip()) #segsites line
        positionsLine = msStream.readline().strip()
        if not positionsLine.startswith("positions:"):
            sys.exit("Malformed ms-style output file. AAAARRRRGGHHH!!!!!\n")
        repStr.append(positionsLine) #positions line

        for i in range(numSamples):
            currLine = msStream.readline()
            repStr.append(currLine.strip())
        line = msStream.readline()
        #advance to the next non-empty line or EOF
        while line and line.strip() == "":
            line = msStream.readline()
        repStr = "\n".join(repStr)
        repLs.append(repStr)
    msStream.close()

    return numSamples, repLs

numSamples1, repLs1 = readAllMSRepsFromFile(msFileName1)
numSamples2, repLs2 = readAllMSRepsFromFile(msFileName2)
if numSamples1 != numSamples2:
    sys.exit("sample size differs between %s (%s) and %s (%s). AAAARRRRGGGGHHHHHHHH!\n" %(msFileName1, numSamples1, msFileName2, numSamples2))
print "./msStyle %s %s\nblah\n" %(numSamples1, numReps1+numReps2)
if shuffle == "shuffle":
    random.shuffle(repLs1)
    random.shuffle(repLs2)
print "\n".join(repLs1[:numReps1])
print "\n".join(repLs2[:numReps2])
