#!/usr/bin/env python
import sys, gzip

#example: python ~/kerncode/msUtils/getMSRepSubrange.py test.msout 1 2 | less
msFileName, firstRepNumber, lastRepNumber = sys.argv[1:] #firstRepNumber and lastRepNumber are zero-based incidces of the first and last reps we want to print
firstRepNumber, lastRepNumber = int(firstRepNumber), int(lastRepNumber)
if lastRepNumber < firstRepNumber:
    sys.exit("lastRepNumber must be >= firstRepNumber. AAAAARRRRRGGGGGGGHHHHHH!!!\n")

if msFileName.endswith(".gz"):
    msStream = gzip.open(msFileName)
elif msFileName == "stdin":
    msStream = sys.stdin
else:
    msStream = open(msFileName)

header = msStream.readline().strip().split()
program,numSamples,numSims = header[:3]
if len(header) > 3:
    otherParams = " " + " ".join(header[3:])
else:
    otherParams = ""
numSamples,numSims = int(numSamples),int(numSims)
sys.stdout.write("./msStyle %s %s%s\nblah" %(numSamples, lastRepNumber-firstRepNumber+1, otherParams))

processedSims = 0
#advance to first simulation
line = msStream.readline()
while not line.strip().startswith("//"):
    line = msStream.readline()
while line:
    if not line.strip().startswith("//"):
        sys.exit("Malformed ms-style output file: read '%s' instead of '//'. AAAARRRRGGHHH!!!!!\n" %(line.strip()))
    repStr = "\n\n//\n"
    repStr += msStream.readline() #segsites line
    positionsLine = msStream.readline()
    if not positionsLine.startswith("positions:"):
        sys.exit("Malformed ms-style output file. AAAARRRRGGHHH!!!!!\n")
    repStr += positionsLine #positions line

    for i in range(numSamples):
        currLine = msStream.readline()
        repStr += currLine
    if processedSims >= firstRepNumber:
        sys.stdout.write(repStr.rstrip())
    if processedSims == lastRepNumber:
        break
    processedSims += 1
    line = msStream.readline()
    #advance to the next non-empty line or EOF
    while line and line.strip() == "":
        #if processedSims > firstRepNumber:
        #    print line.strip()
        line = msStream.readline()

if msFileName != "stdin":
    msStream.close()
