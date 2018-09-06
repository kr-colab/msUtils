import sys, gzip, bisect

def sortedFlankingPositionsByDistToTargSite(targetPos, flankingPositionsToExamine, desiredNumPositions, physLen):
    i=1
    sortedFlankingPositions = []

    while len(sortedFlankingPositions) < desiredNumPositions:
        lPos = targetPos-i
        rPos = targetPos+i
        if lPos >= 0 and lPos in flankingPositionsToExamine:
            sortedFlankingPositions.append(lPos)
        if rPos < physLen and rPos in flankingPositionsToExamine and len(sortedFlankingPositions) < desiredNumPositions:
            sortedFlankingPositions.append(rPos)
        i += 1

    return sortedFlankingPositions

def getNearestEmptyPositions(donorPos, snpCountAtPos, physLen):
    numColliders = snpCountAtPos[donorPos]-1

    freeSlots = []
    lPos = donorPos - 1
    rPos = donorPos + 1
    while len(freeSlots) < numColliders:
        if lPos >= 0:
            if snpCountAtPos[lPos] == 0:
                freeSlots.append(lPos)
            lPos -= 1
        if rPos <= physLen-1:
            if snpCountAtPos[rPos] == 0:
                freeSlots.append(rPos)
            rPos += 1

    return freeSlots

def resolveCollision(donorPos, snpCountAtPos, physLen):
    for recipientPos in getNearestEmptyPositions(donorPos, snpCountAtPos, physLen):
        snpCountAtPos[recipientPos] += 1
        assert snpCountAtPos[recipientPos] == 1
        snpCountAtPos[donorPos] -= 1

def msPositionsToIntegerPositions(positions, physLen):
    assert physLen >= len(positions)

    snpCountAtPos = {}
    for i in range(physLen):
        snpCountAtPos[i] = 0
    for position in positions:
        intPos = int(physLen*position)
        if intPos == physLen:
            intPos = physLen-1
        snpCountAtPos[intPos] += 1

    collisions = {}
    for pos in snpCountAtPos:
        if snpCountAtPos[pos] > 1:
            collisions[pos] = 1

    midPos = physLen/2
    collisionPositions = []
    midHasCollision=0
    if midPos in collisions:
        collisionPositions.append(midPos)
        midHasCollision=1
    collisionPositions += sortedFlankingPositionsByDistToTargSite(midPos, collisions, len(collisions)-midHasCollision, physLen)
    for pos in collisionPositions:
        resolveCollision(pos, snpCountAtPos, physLen)

    assert max(snpCountAtPos.values()) == 1
    newPositions = [x for x in sorted(snpCountAtPos) if snpCountAtPos[x] > 0]
    assert newPositions[0] >= 0 and newPositions[-1] < physLen

    return newPositions

def msRepToHaplotypeArrayIn(samples, positions, totalPhysLen, positionsFirst=True):
    for i in range(len(samples)):
        assert len(samples[i]) == len(positions)

    positions = msPositionsToIntegerPositions(positions, totalPhysLen)

    hapArrayIn = []
    if positionsFirst:
        for j in range(len(positions)):
            hapArrayIn.append([])
            for i in range(len(samples)):
                hapArrayIn[j].append(samples[i][j])
    else:
        for i in range(len(samples)):
            hapArrayIn.append([])
            for j in range(len(positions)):
                hapArrayIn[i].append(samples[i][j])

    return hapArrayIn, positions

def msOutToHaplotypeArrayIn(msOutputFileName, totalPhysLen, positionsFirst=True):
    if msOutputFileName == "stdin":
        isFile = False
        msStream = sys.stdin
    else:
        isFile = True
        if msOutputFileName.endswith(".gz"):
            msStream = gzip.open(msOutputFileName)
        else:
            msStream = open(msOutputFileName)

    header = msStream.readline()
    program,numSamples,numSims = header.strip().split()[:3]
    numSamples,numSims = int(numSamples),int(numSims)

    hapArraysIn = []
    positionArrays = []
    #advance to first simulation
    line = msStream.readline()
    while not line.strip().startswith("//"):
        line = msStream.readline()
    while line:
        if not line.strip().startswith("//"):
            sys.exit("Malformed ms-style output file: read '%s' instead of '//'. AAAARRRRGGHHH!!!!!\n" %(line.strip()))
        segsitesBlah,segsites = msStream.readline().strip().split()
        segsites = int(segsites)
        if segsitesBlah != "segsites:":
            sys.exit("Malformed ms-style output file. AAAARRRRGGHHH!!!!!\n")

        positionsLine = msStream.readline().strip().split()
        if not positionsLine[0] == "positions:":
            sys.exit("Malformed ms-style output file. AAAARRRRGGHHH!!!!!\n")
        positions = [float(x) for x in positionsLine[1:]]

        samples = []
        for i in range(numSamples):
            sampleLine = msStream.readline().strip()
            if len(sampleLine) != segsites:
                sys.exit("Malformed ms-style output file %s segsites but %s columns in line: %s; line %s of %s samples AAAARRRRGGHHH!!!!!\n" %(segsites,len(sampleLine),sampleLine,i,numSamples))
            samples.append(sampleLine)
        if len(samples) != numSamples:
            raise Exception
        hapArrayIn, positions = msRepToHaplotypeArrayIn(samples, positions, totalPhysLen, positionsFirst=positionsFirst)
        hapArraysIn.append(hapArrayIn)
        positionArrays.append(positions)
        line = msStream.readline()
        #advance to the next non-empty line or EOF
        while line and line.strip() == "":
            line = msStream.readline()
        #sys.stderr.write("finished rep %d\n" %(len(hapArraysIn)))
    if len(hapArraysIn) != numSims:
        sys.exit("Malformed ms-style output file: %s of %s sims processed. AAAARRRRGGHHH!!!!!\n" %(len(hapArraysIn), numSims))

    if isFile:
        msStream.close()
    return hapArraysIn, positionArrays

def msOutToHaplotypeMatrices(msOutputFileName, totalPhysLen):
    return msOutToHaplotypeArrayIn(msOutputFileName, totalPhysLen, positionsFirst=False)

def windowHaps(hapArraySamplesFirst, positionArray, winStart, winEnd):
    indicesToKeep = [i for i in range(len(positionArray)) if positionArray[i] >= winStart and positionArray[i] <= winEnd]
    windowHapArray = []
    for sample in hapArraySamplesFirst:
        windowHapArray.append([])
        for i in indicesToKeep:
            windowHapArray[-1].append(sample[i])
    windowPositions = [positionArray[i] for i in indicesToKeep]
    return windowHapArray, windowPositions

def msWinStr(hapArraysSamplesFirst, positionArrays, winStart, winEnd):
    outStr = "./windowedMSOutput %s %s\nblah\n" %(len(hapArraysSamplesFirst[0]), len(hapArraysSamplesFirst))
    for i in range(len(hapArraysSamplesFirst)):
        currHapArray, currPositions = windowHaps(hapArraysSamplesFirst[i], positionArrays[i], winStart, winEnd)
        currPositions = [(pos-winStart)/(winEnd-winStart+1.0) for pos in currPositions]
        outStr += "\n//\nsegsites: %s\n" %(len(currPositions))
        outStr += "positions: " + " ".join([str(x) for x in currPositions]) + "\n"
        for sample in currHapArray:
            outStr += "".join(sample) + "\n"
    return outStr
