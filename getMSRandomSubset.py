#!/usr/bin/env python
import sys, gzip, random

msFileName, numReps, header, headerNumReps = sys.argv[1:]
if headerNumReps != "actual":
    if headerNumReps < 0:
        sys.exit("headerNumReps must be < 0 or 'na' or 'actual'. AAAARRRGGGHHH!!!\n")
    elif headerNumReps == "na":
        assert header == "no_header"
    else:
        headerNumReps = int(headerNumReps)
numReps = int(numReps)
if not header in ["header", "no_header"]:
    sys.exit("header must be set to either 'header' or 'no_header'. AAAARRRRGGGGHHHHHHHHHH!!!!\n")

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

numSamples, repLs = readAllMSRepsFromFile(msFileName)
if header == "header":
    if headerNumReps == "actual":
        print "./msStyle %s %s\nblah" %(numSamples, numReps)
    else:
        print "./msStyle %s %s\nblah" %(numSamples, headerNumReps)

if numReps > 0:
    random.shuffle(repLs)
    print "\n".join(repLs[:numReps])
