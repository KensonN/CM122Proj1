import sys
from copy import deepcopy

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
    possibleIndices = dict()
    substitutions = set()

    def __init__(self, genomeFile, readFile):
        self.genomeFile = genomeFile
        self.readFile = readFile
        self.importGenome()
        self.importReads()
        self.mapGenome()
        self.writeGenomeMap()
        for i in range(0, self.genomeSize):
            self.coverageMap[i] = 0
            self.possibleIndices[i] = []
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
                mutation.append((i+alignment, readString[i], genomeString[i]))
        return mutation

    def align(self):
        for read in self.reads:
            # print(i[0], i[1], "READ")
            positions = [[-10]] * (len(read[1])-self.kmerSize+1)
            # print(i[0], i[1], len(i[1])-self.kmerSize+1)
            for k in range(0, len(read[1])-self.kmerSize+1):
                keys = list(self.genomeMap.keys())
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
            last = None
            # print(positions)
            finalPositions = []
            # for i in range(len(positions)):
            #     if positions[i] == [-10]:
            #         print(i, end=", ")
            # print()

            counter = 0
            for k in positions:
                for j in k:
                    # print("j=", j, end="| ")
                    if last == None:
                        last = j
                        finalPositions.append(j)
                        counter = 0
                        # print(j, "0")
                    elif abs((j-last)) <= 2 and j != -10: 
                        last = j
                        finalPositions.append(j)
                        counter += 1
                        # print(j, "1")
                    elif counter > 8 and j == -10:
                        last += 1
                        finalPositions.append(last)
                        counter+=1
            # print(finalPositions)
            first = None
            finalPositionsTwo = deepcopy(finalPositions)
            modifiedRead = read[1]
            for k in finalPositions:
                if first == None:
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
                        if deletionLoc in self.deletions:
                            self.deletions[deletionLoc] += 1
                        else:
                            self.deletions[deletionLoc] = 1
                        # print(int((i+first)/2))
                if k == first:
                    # print("insertion: ")
                    insertionLoc = first-1
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

                    if insertionLoc in self.insertions:
                        self.insertions[insertionLoc] += 1
                    else:
                        self.insertions[insertionLoc] = 1
                    first = k
                first = k
            # print(finalPositionsTwo)
            # print(read[1])
            # print(self.genome[262:262+50])
            if [-10] != finalPositionsTwo and self.hammingDistance(modifiedRead, self.genome[finalPositionsTwo[0]:finalPositionsTwo[0]+len(modifiedRead)]) <= self.allowedErrors:
                
                mutations = (self.getSubstitutions(modifiedRead, self.genome[finalPositionsTwo[0]:finalPositionsTwo[0]+len(modifiedRead)], finalPositionsTwo[0]))
                # if len(mutations) > 0:
                #     print(read[1], "READ")
                #     print(modifiedRead, "MODIFIED")
                #     print(self.genome[finalPositionsTwo[0]:finalPositionsTwo[0]+len(modifiedRead)], "GENOME")
                for i in mutations:
                    self.substitutions.add(i)
                    print(i)
                for position in range(len(finalPositionsTwo)):
                    for i in range(self.kmerSize):
                        # print(finalPositionsTwo[position] + i, modifiedRead[i+position], self.genome[finalPositionsTwo[position] + i]  ,end=" ")
                        # print()
                        self.possibleIndices[finalPositionsTwo[position] + i].append(modifiedRead[i+position])
                        pass
                    # print()
                    # print("---")
                    pass
        # print(self.possibleIndices)
        # print(found, total)
        # print(self.coverageMap)
        print("deletions: ", self.deletions)
        print("insertions: ", self.insertions)
        print("subsitutions:", self.substitutions)

referenceGenome = "sample_1000\sample_1000_reference_genome.fasta"
# referenceGenome = "project1a_10000_reference_genome.fasta"

reads = "sample_1000\sample_1000_no_error_single_reads.fasta"
# reads = "sample_1000\sample_1000_with_error_single_reads.fasta"
# reads = "sample_1000\sample_1000_with_error_paired_reads.fasta"
# reads = "sample_1000\DELETION.fasta"
# reads = "sample_1000\INSERTION.fasta"
# reads = "sample_1000\REDUCED.fasta"
# reads = "project1a_10000_with_error_paired_reads.fasta"

aligner = readAligner(genomeFile=referenceGenome, readFile=reads)
aligner.align()