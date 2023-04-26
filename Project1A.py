import sys
from copy import deepcopy
import time

class readAligner():
    kmerSize = 10
    genomeFile = None
    readFile = None
    reads = None
    genome = ""
    genomeMap = None
    allowedErrors = 2
    insertions = dict()
    deletions = dict()
    coverageMap = dict()
    # substitutions = set()
    positionMap = dict()
    subsitutions = dict()
    positions = []

    def __init__(self, genomeFile, readFile):
        self.genomeFile = genomeFile
        self.readFile = readFile
        self.importGenome()
        self.importReads()
        self.mapGenome()
        self.writeGenomeMap()
        for i in range(0, self.genomeSize):
            self.coverageMap[i] = 0
            self.positionMap[i] = []
        # self.coverageMap[-10] = 0

    def importGenome(self):
        with open(self.genomeFile) as f:
            lines = f.read().splitlines()
            # for i in range(len(lines)):
            #     lines[i] = lines[i].strip()
            info = lines[0]
            genome = lines[1:]
            fullGenome = ""
            for i in genome:
                fullGenome += i
            for i in fullGenome:
                if i == '\n':
                    # print(i)
                    pass
            self.genome = fullGenome
            self.genomeSize = len(fullGenome)
        pass

    def importReads(self):
        with open(self.readFile) as f:
            lines = f.read().splitlines()
            # print(lines[0][6:])
            reads = [None]
            readType = ""
            if "/" in lines[0][6:]:
                readType = "paired"
            else:
                readType = "single"
            for i in range(0, len(lines), 2):
                reads.append([lines[i], lines[i+1]])
                pass
            reads = reads[1:]
            self.reads = reads
            # print(self.reads)
        pass

    def mapGenome(self):
        genomeMap = dict()
        for i in range(0, len(self.genome)-self.kmerSize+1):
            if self.genome[i:i+self.kmerSize] in genomeMap:
                genomeMap[self.genome[i:i+self.kmerSize]].append(i)
            else:
                genomeMap[self.genome[i:i+self.kmerSize]] = [i]
        pass
        self.genomeMap = deepcopy(genomeMap)

    def writeGenomeMap(self):
        f = open("genomeMap.txt", "w")
        for i in self.genomeMap:
            # print(i, self.genomeMap[i])
            f.write(i)
            f.write(":")
            for j in self.genomeMap[i]:
                f.write(str(j))
                f.write(",")
            f.write('\n')


    def hammingDistance(self, str1, str2):
        hamming = 0
        if len(str1) != len(str2):
            print(str1, "s1")
            print(str2, "s2")
            raise ValueError("Input strings must have the same length")
        for i in range(len(str1)):
            if str1[i] != str2[i]:
                hamming += 1
        return hamming

    # def hamming_distance(self, s1, s2):
    #     if len(s1) != len(s2):
    #         raise ValueError("Input strings must have the same length")
    #     distance = 0
    #     for i in range(len(s1)):
    #         if s1[i] != s2[i]:
    #             distance += 1
    #     return distance

    def getSubstitutions(self, genomeString, readString, alignment):
        if len(genomeString) != len(readString):
            print(genomeString)
            print(readString)
            raise ValueError("Input strings must have the same length")
        mutation = []
        for i in range(len(genomeString)):
            if genomeString[i] != readString[i] and genomeString[i] != ' ':
                if (i+alignment, readString[i], genomeString[i]) not in self.subsitutions:
                    self.subsitutions[(i+alignment, readString[i], genomeString[i])] = 1
                else:
                    self.subsitutions[(i+alignment, readString[i], genomeString[i])] +=1
                mutation.append((i+alignment, readString[i], genomeString[i]))
        return mutation

    def align(self):
        firstPair = True
        lastReadPos = 0
        startTime = time.perf_counter()
        lastCurrentTime = 0
        for read in self.reads:
            # print(lastReadPos)
            currentTime = time.perf_counter() - startTime
            print(read[0], read[1], "READ", lastReadPos, currentTime, currentTime - lastCurrentTime)
            lastCurrentTime = currentTime
            positions = [[-10]] * (len(read[1])-self.kmerSize+1)
            # print(i[0], i[1], len(i[1])-self.kmerSize+1)
            for k in range(0, len(read[1])-self.kmerSize+1):
                keys = list(self.genomeMap.keys())[lastReadPos:]
                for j in keys:
                    # print(j)
                    if self.hammingDistance(read[1][k:k+self.kmerSize], j) == 0:
                        if positions[k] == [-10]:
                            # print("if")
                            positions[k] = deepcopy(self.genomeMap[j])
                        else:
                            # print("else")
                            positions[k] += deepcopy(self.genomeMap[j])
                # total += 1
            # print(positions, "pos")
            finalPositions = []
            # for i in range(len(positions)):
            #     if positions[i] == [-10]:
            #         print(i, end=", ")
            # print()
            last = None
            counter = 0
            for k in positions:
                for j in k:
                    # print("j=", j, end="| ")
                    if last == None:
                        if j == -10:
                            finalPositions.append(j)
                        else:
                            last = j
                            finalPositions.append(j)
                        # print(j, "0")
                    elif j!= -10:
                        if abs((j-last)) <= 2:
                            last = j
                            finalPositions.append(j)
                            counter = 0
                    # elif abs((j-last)) <= 2 and j != -10: 
                    #     last = j
                    #     finalPositions.append(j)
                    #     counter = 0
                    #     # print(j, "1")
                    elif j == -10:
                        if counter <= self.kmerSize:
                            last += 1
                            finalPositions.append(last)
                            counter+=1
                    # elif counter <= self.kmerSize and j == -10:
                    #     last += 1
                    #     finalPositions.append(last)
                    #     counter+=1
            last = None
            counter = 0
            # for k in range(1, len(positions)):
            #     for j in k:
            #         # print("j=", j, end="| ")
            #         if last == None:
            #             if j == -10:
            #                 finalPositions.append(j)
            #             else:
            #                 last = j
            #                 finalPositions.append(j)
            #             # print(j, "0")
            #         elif abs((j-last)) <= 2 and j != -10: 
            #             last = j
            #             finalPositions.append(j)
            #             counter = 0
            #             # print(j, "1")
            #         elif counter <= self.kmerSize and j == -10:
            #             last += 1
            #             finalPositions.append(last)
            #             counter+=1
                # print(counter, j)
            # print(finalPositions, "finalPos")
            # print("---")
            first = None
            finalPositionsTwo = deepcopy(finalPositions)
            modifiedRead = read[1]
            for k in finalPositions:
                if first == None:
                    if k != -10:
                        first = k
                    continue
                if k != first + 1:
                    diff = k - first
                    if diff == 2:
                        deletionLoc = int((k+first)/2)
                        self.coverageMap[deletionLoc] += 1
                        # finalPositionsTwo.insert(finalPositionsTwo.index(first)+1, first+1)
                        delIndex = finalPositions.index(first+2)
                        # print(read[1][:delIndex])
                        # print(self.genome[finalPositions[0] + delIndex: finalPositions[0] + delIndex + 10])
                        modifiedRead = modifiedRead[:delIndex] + self.genome[finalPositions[0] + delIndex] + modifiedRead[delIndex:]
                        # print(self.genome[finalPositions[0]:finalPositions[0]+50], "genome")
                        # print(read[1], "read")
                        # print("---")
                        # print(self.genome[finalPositions[0]:finalPositions[0]+50], "genome")
                        # print(modifiedRead, "modified read")
                        deletion = (deletionLoc, self.genome[finalPositions[0] + delIndex])
                        if deletion in self.deletions:
                            self.deletions[deletion] += 1
                        else:
                            self.deletions[deletion] = 1
                        # print(int((i+first)/2))
                if k == first:
                    # print("insertion: ")
                    insertionLoc = first-1
                    # print(insertionLoc)
                    # self.coverageMap[insertionLoc+1] -= 1
                    finalPositionsTwo.remove(insertionLoc+1)
                    insertionIndex = finalPositions.index(insertionLoc+1)
                    modifiedRead = modifiedRead[:insertionIndex] + modifiedRead[insertionIndex+1:]

                    # print(self.genome[finalPositions[0]:finalPositions[0]+50])
                    # print(modifiedRead)
                    # print(read[1])
                    # print("---")
                    # print(self.genome[finalPositions[0]:finalPositions[0]+50])
                    # print(modifiedRead)
                    # print(modifiedRead[:insertionIndex])
                    # print(read[1][insertionIndex])
                    # print(modifiedRead[insertionIndex+1:])   
                    insertion = (insertionLoc, read[1][insertionIndex])
                    if insertion in self.insertions:
                        self.insertions[insertion] += 1
                    else:
                        self.insertions[insertion] = 1
                    first = k
                first = k
            self.positions.append(finalPositionsTwo)


            lastPos = None
            if len(finalPositionsTwo) > 0:
                for pos in finalPositionsTwo:
                    if pos != -10:
                        self.coverageMap[pos] += 1
                        lastPos = pos
                        # print(pos)
                if lastPos != None:
                    for j in range(1, self.kmerSize):
                        self.coverageMap[lastPos+j] += 1
                    
                    if firstPair == True:
                        lastReadPos = lastPos+self.kmerSize-1


                stringAlignment = 0
                genomeString = ""
                print(lastReadPos, lastReadPos)
                print(positions)
                print(finalPositions)
                print(finalPositionsTwo)
                if finalPositionsTwo[0] == -10 and finalPositionsTwo[-1] == -10:
                    stringAlignment = 999
                elif finalPositionsTwo[0] == -10:
                    genomeString = self.genome[finalPositionsTwo[-1]-len(modifiedRead)+self.kmerSize:finalPositionsTwo[-1]+self.kmerSize]
                    stringAlignment = self.hammingDistance(modifiedRead, self.genome[finalPositionsTwo[-1]-len(modifiedRead):finalPositionsTwo[-1]])
                else:
                    genomeString = self.genome[finalPositionsTwo[0]:finalPositionsTwo[0]+len(modifiedRead)]
                    stringAlignment = self.hammingDistance(modifiedRead, genomeString)

                stringAlignment = self.hammingDistance(modifiedRead, genomeString)
                
                if [-10] != finalPositionsTwo and stringAlignment <= self.allowedErrors:
                    
                    mutation = self.getSubstitutions(modifiedRead, self.genome[finalPositionsTwo[0]:finalPositionsTwo[0]+len(modifiedRead)], finalPositionsTwo[0])

            if firstPair == True:
                # print("T")
                firstPair = False
            elif firstPair == False:
                # print("F")
                firstPair = True
                lastReadPos = 0


        print("test")

        print("Deletions:")
        for i in self.deletions:
            print(i, self.deletions[i], self.coverageMap[i[0]])
        print("---")
        # print("deletions: ", self.deletions)
        print("Insertions: ", self.insertions)
        for i in self.insertions:
            print(i, self.insertions[i], self.coverageMap[i[0]])
        print("---")
        print("Substitutions:")
        print(self.subsitutions)
        for i in self.subsitutions:
            print(i, self.subsitutions[i], self.coverageMap[i[0]])
        print("---")
        # print(self.positions)
        # for i in range(0, len(self.positions), 2):
        #     print(self.positions[i])
        #     print(self.positions[i+1])
        #     print("---")
        # print(len(self.positions))
        print("Deletions:", self.deletions)
        for i in self.deletions:
            if self.deletions[i] / self.coverageMap[i[0]] > 0.25 and self.coverageMap[i[0]] > 1 and self.deletions[i] > 1:
                print(i, self.deletions[i], self.coverageMap[i[0]])
        print("---")
        # print("deletions: ", self.deletions)
        print("Insertions: ", self.insertions)
        for i in self.insertions:
            if self.insertions[i] / self.coverageMap[i[0]] > 0.25 and self.coverageMap[i[0]] > 1 and self.insertions[i] > 1:
                print(i, self.insertions[i], self.coverageMap[i[0]])
        print("---")
        print("Substitutions:")
        # print(self.subsitutions)
        for i in self.subsitutions:
            if self.subsitutions[i] / self.coverageMap[i[0]] > 0.25:
                print(i, self.subsitutions[i], self.coverageMap[i[0]])
        # print(self.coverageMap)
        # print("Subsitutions:", sorted(self.substitutions))

referenceGenome = "sample_1000\sample_1000_reference_genome.fasta"
# referenceGenome = "project1a_10000_reference_genome.fasta"

# reads = "sample_1000\sample_1000_no_error_single_reads.fasta"
# reads = "sample_1000\sample_1000_with_error_single_reads.fasta"
# reads = "sample_1000\sample_1000_with_error_paired_reads.fasta"
# reads = "sample_1000\sample_1000_no_error_paired_reads REDUCED.fasta"
# reads = "sample_1000\sample_1000_no_error_paired_reads.fasta"
reads = "sample_1000\DELETION.fasta"
# reads = "sample_1000\INSERTION.fasta"
# reads = "sample_1000\REDUCED.fasta"
# reads = "project1a_10000_with_error_paired_reads.fasta"
# reads = "project1a_10000_with_error_paired_reads_TEST.fasta"

aligner = readAligner(genomeFile=referenceGenome, readFile=reads)
aligner.align()