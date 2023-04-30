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
    # positionMap = dict()
    substitutions = dict()
    finalInsertions = dict()
    finalDeletions = dict()
    finalSubstitutions = dict()
    positions = []
    stepSize = 1

    def __init__(self, genomeFile, readFile):
        self.genomeFile = genomeFile
        self.readFile = readFile
        self.importGenome()
        self.importReads()
        self.mapGenome()
        self.writeGenomeMap()
        for i in range(0, self.genomeSize):
            self.coverageMap[i] = 0
            # self.positionMap[i] = []
        # self.coverageMap[-10] = 0
        self.align()
        self.writeResults()

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
                if (i+alignment, readString[i], genomeString[i]) not in self.substitutions:
                    self.substitutions[(i+alignment, readString[i], genomeString[i])] = 1
                else:
                    self.substitutions[(i+alignment, readString[i], genomeString[i])] +=1
                mutation.append((i+alignment, readString[i], genomeString[i]))
        return mutation

    def getStreak(self, positions):
        avg = 0
        counter = 0
        for i in positions:
            if i != [-10] and len(i) == 1:
                avg += i[0]
                counter+=1
        avg /= counter
        positionsReduced = []
        for i in positions:
            if i!= [-10] and len(i) > 1:
                index = 0
                diff = 999999
                for count, arrayItem in enumerate(i):
                    if abs(arrayItem - avg) < diff:
                        diff = abs(arrayItem - avg)
                        index = count
                positionsReduced.append(i[index])
            else:
                positionsReduced.append(i[0])

        streak = False
        currentStreakSum= 0
        currentStreakStartIndex = 0
        currentStreakCount = 0
        streakAvgs = []
        streakIndexes = []
        streakCounts = []
        last = None
        for count, i in enumerate(positionsReduced):
            # print("count:", count)
            # print(last)
            if i != -10:
                if streak == False:
                    # print(i)
                    streak = True
                    currentStreakStartIndex = count
                    currentStreakSum = i
                    currentStreakCount = 1
                elif streak == True:
                    if abs(i - last)/self.stepSize < 2:
                        currentStreakSum += i
                        currentStreakCount += 1
                        if count == (len(positionsReduced)-1):
                            streak = False
                            streakAvgs.append(currentStreakSum/currentStreakCount)
                            streakIndexes.append(currentStreakStartIndex)
                            streakCounts.append(currentStreakCount)
                    else:
                        # print("else")
                        streakAvgs.append(currentStreakSum/currentStreakCount)
                        streakIndexes.append(currentStreakStartIndex)
                        streakCounts.append(currentStreakCount)
                        currentStreakStartIndex = count
                        currentStreakSum = i
                        currentStreakCount = 1
                last = i
            elif i == -10:
                if streak == True:
                    streak = False
                    streakAvgs.append(currentStreakSum/currentStreakCount)
                    streakIndexes.append(currentStreakStartIndex)
                    streakCounts.append(currentStreakCount)



        maxStreakCountIndex = streakCounts.index(max(streakCounts))
        maxStreakStartIndex = streakIndexes[maxStreakCountIndex]

        # print(streakCounts)
        # print(streakIndexes)
        # print(maxStreakStartIndex)

        streakIndexes.pop(maxStreakCountIndex)
        streakCounts.pop(maxStreakCountIndex)
        # print(positionsReduced, "reduced")
        for count, i in enumerate(streakIndexes):
            # print(i, i+streakCounts[count])
            for j in range(i, i+streakCounts[count]):
                # print(j, end=" ")
                # print()
                # print(abs(abs(maxStreakStartIndex - j) - abs(positionsReduced[maxStreakStartIndex] - positionsReduced[j])))

                # print(abs(positionsReduced[maxStreakStartIndex] - positionsReduced[j]))
                # print(abs(maxStreakStartIndex - j) - abs(positionsReduced[maxStreakStartIndex] - positionsReduced[j]))
                # print(positionsReduced[maxStreakStartIndex])
                # print(j)
                # print("---")

                if abs(abs(positionsReduced[maxStreakStartIndex] - positionsReduced[j]) - abs(maxStreakStartIndex - j)*self.stepSize) > 3:

                # if (abs(abs(maxStreakStartIndex - j) - abs(positionsReduced[maxStreakStartIndex] - positionsReduced[j])))/self.stepSize > 3:
                    positionsReduced[j] = -10
                    # print(abs(maxStreakStartIndex - j))
                    # print(abs(positionsReduced[maxStreakStartIndex] - positionsReduced[j]))
                    # print(abs(maxStreakStartIndex - j) - abs(positionsReduced[maxStreakStartIndex] - positionsReduced[j]))
                    # print("---")
                pass
                # print("--")
            # print("---")
            pass        
        pass
        # print('\n\n\n')
        # print(positionsReduced, "posRedux")
        return positionsReduced
        returnArray = []
        for i in positionsReduced:
            returnArray.append([i])
        # print(returnArray)
        return returnArray

    def align(self):
        # firstPair = True
        lastReadPos = 0
        startTime = time.perf_counter()
        lastCurrentTime = 0
        for count, read in enumerate(self.reads):
            # print(lastReadPos)
            # print(read)
            readSplit = read[0].split("/")
            # lastReadPos = 0
            if len(readSplit) > 1:
                if read[0].split("/")[1] == '1':
                    lastReadPos = 0
            else:
                lastReadPos = 0

            currentTime = time.perf_counter() - startTime
            # print(read[0], '\n', "readPos:", lastReadPos, '\n', "currentTime:", round(currentTime, 3), '\n', "readTime:",  round(currentTime - lastCurrentTime, 3))

            # if count == len(self.reads) - 1:
            #     print(read[0], lastReadPos, round(currentTime, 3), round(currentTime - lastCurrentTime, 3), end='\n')
            # else:
            #     print(read[0], lastReadPos, round(currentTime, 3), round(currentTime - lastCurrentTime, 3), end='\r')

            insertionState = False
            deletionState = False

            lastCurrentTime = currentTime
            positions = [[-10]] * (int((len(read[1])/1-self.kmerSize)/self.stepSize) + 1)
            for kmerPos in range(0, len(read[1])-self.kmerSize+1, self.stepSize):
                keys = list(self.genomeMap.keys())[lastReadPos:]
                counter = 0
                for key in keys:
                    if read[1][kmerPos:kmerPos+self.kmerSize] == key:
                        if positions[int(kmerPos/self.stepSize)] == [-10]:
                            positions[int(kmerPos/self.stepSize)] = deepcopy(self.genomeMap[key])
                        else:
                            positions[int(kmerPos/self.stepSize)] += deepcopy(self.genomeMap[key])

            if len(positions) == 0:
                continue

            positionsFiltered = self.getStreak(positions)
            
            finalPositions = []
            last = None
            counter = 0
            streak = 0
            for pos, posItem in enumerate(positionsFiltered):
                if last == None:
                    if posItem == -10:
                        finalPositions.append(posItem)
                        counter += 1
                    else:
                        last = posItem
                        finalPositions.append(posItem)
                        streak+=1
                        if pos > 0:
                            for l in range(1, counter+1):
                                finalPositions[pos-l] = posItem - l*self.stepSize
                                streak+=1
                                pass
                elif posItem != -10:
                    if abs(posItem-last) <= 5:
                        last = posItem
                        finalPositions.append(posItem)
                        streak += 1
                        counter = 0
                    else:
                        if streak < 2:
                            for l in range(0, counter+1):
                                # print(pos-l-1, "posl")
                                finalPositions[pos-l-1] = -10
                                pass
                            last = posItem
                            finalPositions.append(posItem)
                            counter = 0
                        finalPositions.append(last+1*self.stepSize)
                        last +=1*self.stepSize
                elif posItem == -10:
                    last += 1*self.stepSize
                    finalPositions.append(last)
                    counter+=1
            last = None
            counter = 0

            #fill in -10s            
            if finalPositions[0] == -10:
                for count, item in enumerate(finalPositions):
                    if item != -10:
                        for i in range(1, count+1):
                            finalPositions[count - i] = item - i
                        break
            
            # print(finalPositions, "finalPositions")

            first = None
            finalPositionsTwo = deepcopy(finalPositions)
            modifiedRead = read[1]

            if -10 in finalPositions:
                lastReadPos = 0
                print()
                print("skipping", read)
                print()
                continue
                # f = open("skips.txt", "w")
                # f.write("Skipping: ")
                # f.write(str(read[0]))
                # f.write(" ")
                # f.write(str(read[1]))
                # f.write('\n')
                # continue

            print(positions, "pos")
            print(positionsFiltered, "posFilt")
            print(finalPositions, "final")
            print("---")

            
            for k in finalPositions:
                if first == None:
                    if k != -10:
                        first = k
                    continue
                # diff = (k - first)
                # print("diff:", diff)
                if k != first + 1*self.stepSize:
                    diff = (k - first)
                    print("diff:", diff)
                    # print(finalPositions)
                    
                    # print(positions)
                    # print(positionsFiltered)
                    # print(finalPositions)
                    if diff == 1+self.stepSize:
                        # print(read)
                        # print(finalPositions)
                        # deletionLoc = int((k+first)/2)
                        # self.coverageMap[deletionLoc] += 1
                        # delIndex = finalPositions.index(first+2)
                        # print(delIndex)
                        # modifiedRead = modifiedRead[:delIndex] + self.genome[finalPositions[0] + delIndex] + modifiedRead[delIndex:]

                        # print(read[0])
                        print("k", k, "first", first)
                        # deletionLoc = int((k+first)/2)
                        deletionLoc = k - self.stepSize
                        print("deletionLoc:", deletionLoc)
                        # self.coverageMap[deletionLoc] += 1
                        delIndex = (finalPositions.index(first)+1)*self.stepSize
                        print(delIndex, "delIndex")
                        print(modifiedRead[:delIndex])
                        modifiedRead = modifiedRead[:delIndex] + self.genome[deletionLoc] + modifiedRead[delIndex:]


                        if deletionLoc < 0:
                            print(read)
                            print(positions, "pos")
                            print(finalPositions, "final")
                            print(finalPositionsTwo, "finalTwo")
                            print("---")
                            break
                        deletion = (deletionLoc, self.genome[finalPositions[0] + delIndex])
                        if deletion in self.deletions:
                            self.deletions[deletion] += 1
                        else:
                            self.deletions[deletion] = 1
                        deletionState = True
                    elif diff == self.stepSize-1:
                        # print(positions)
                        # print(positionsFiltered)
                        # print(finalPositions)

                        insertionLoc = first-1

                        # print("insertionLoc", insertionLoc)

                        if insertionLoc < 0:
                            print(read)
                            print(positions, "pos")
                            print(finalPositions, "final")
                            print(finalPositionsTwo, "finalTwo")
                            print("---")
                            break
                    
                    # print(insertionLoc)
                    # print(read)
                    # print(positions, "pos")
                    # print(finalPositions, "final")
                    # print(finalPositionsTwo, "finalTwo")
                    # print("---")
                    # self.coverageMap[insertionLoc] -= 1

                        # print("is in list", insertionLoc+1 in finalPositionsTwo)
                        finalPositionsTwo.remove(insertionLoc+1)
                        insertionIndex = finalPositions.index(insertionLoc+1)
                        modifiedRead = modifiedRead[:insertionIndex] + modifiedRead[insertionIndex+1:]

                        # print("insertion", insertionLoc)
                        # print(read[1], "read")
                        # print(self.genome[finalPositions[0]:finalPositions[0]+50], "genome")
                        # print(modifiedRead ,"modified read")

                        insertion = (insertionLoc, read[1][insertionIndex])
                        if insertion in self.insertions:
                            self.insertions[insertion] += 1
                        else:
                            self.insertions[insertion] = 1
                        first = k
                        insertionState = True

                first = k

            # lastPos = None
            if len(finalPositionsTwo) > 0:
            #     for pos in finalPositionsTwo:
            #         for j in range(self.stepSize):
            #             self.coverageMap[j+pos] += 1
            #             lastPos = pos

            #     if lastPos != None:
            #         for j in range(1, self.kmerSize):
            #             # print(lastPos+j)
            #             self.coverageMap[lastPos+j] += 1
            #             pass
                for pos in range(finalPositionsTwo[0], finalPositionsTwo[-1]+self.kmerSize):
                    self.coverageMap[pos] += 1
                
                # if 413 in finalPositionsTwo:
                #     for pos in finalPositionsTwo:
                #         if pos != -10:
                #             # print(pos)
                #             lastPos = pos
                    # if lastPos != None:
                    #     for j in range(1, self.kmerSize):
                    #         # print(lastPos+j)
                    #         # self.coverageMap[lastPos+j] += 1
                    #         # pass
                    #         # print(lastPos+j)
                    #         pass
                    #     pass
                    # pass
                    # print(self.coverageMap[414])
                    

                lastReadPos = finalPositionsTwo[0]
                stringAlignment = 0
                genomeString = ""

                if finalPositionsTwo[0] == -10 and finalPositionsTwo[-1] == -10:
                    stringAlignment = 999
                elif finalPositionsTwo[0] == -10:
                    genomeString = self.genome[finalPositionsTwo[-1]-len(modifiedRead)+self.kmerSize:finalPositionsTwo[-1]+self.kmerSize]
                    stringAlignment = self.hammingDistance(modifiedRead, genomeString)
                else:
                    genomeString = self.genome[finalPositionsTwo[0]:finalPositionsTwo[0]+len(modifiedRead)]

                    print(read[1], "read")
                    print(genomeString, "genome")
                    print(modifiedRead ,"modified read")
                    print(stringAlignment, "stralign")
                    print(self.hammingDistance(modifiedRead, genomeString))
                    print("insertion", insertionState)
                    print("deletion", deletionState)

                    stringAlignment = self.hammingDistance(modifiedRead, genomeString)

                # stringAlignment = self.hammingDistance(modifiedRead, genomeString)
                self.allowedErrors=3
                if [-10] != finalPositionsTwo and stringAlignment <= self.allowedErrors:
                    
                    mutation = self.getSubstitutions(modifiedRead, genomeString, finalPositionsTwo[0])
            # print(self.coverageMap[414], read)  
            # if firstPair == True:
            #     # print("T")
            #     firstPair = False
            # elif firstPair == False:
            #     # print("F")
            #     firstPair = True
            #     lastReadPos = 0
        
        print("Deletions:" )
        for i in self.deletions:
            print(i, self.deletions[i], self.coverageMap[i[0]])
        # print("---")
        # print("deletions: ", self.deletions)
        print("Insertions: ")
        for i in self.insertions:
            print(i, self.insertions[i], self.coverageMap[i[0]])
        # print("---")
        # print("Substitutions:")
        print("Subsitutions")
        for i in self.substitutions:
            print(i, self.substitutions[i], self.coverageMap[i[0]])

        print('\n', "---", '\n')

        # print("Deletions:" )
        # for i in self.deletions:
        #     # if self.deletions[i] / self.coverageMap[i[0]] > 0.20 and self.coverageMap[i[0]] > 1 and self.deletions[i] > 1:
        #     if self.coverageMap[i[0]] > 1 and self.deletions[i] > 1:
        #         self.finalDeletions[i] = self.deletions[i]
        #         print(i, self.deletions[i], self.coverageMap[i[0]])
        # # print("---")
        # # print("deletions: ", self.deletions)
        # print("Insertions: ")
        # for i in self.insertions:
        #     if self.coverageMap[i[0]] > 1 and self.insertions[i] > 1:
        #     # if self.insertions[i] / self.coverageMap[i[0]] > 0.20 and self.coverageMap[i[0]] > 1 and self.insertions[i] > 1:
        #         self.finalInsertions[i] = self.insertions[i]
        #         print(i, self.insertions[i], self.coverageMap[i[0]])
        # # print("---")
        # # print("Substitutions:")
        # print("Substitutions")
        # for i in self.substitutions:
        #     if self.substitutions[i] / self.coverageMap[i[0]] > 0.5 and self.coverageMap[i[0]] > 1 and self.substitutions[i] > 1:
        #         self.finalSubstitutions[i] = self.substitutions[i]
        #         print(i, self.substitutions[i], self.coverageMap[i[0]])

    def writeResults(self):
        f = open("mutations.txt", "w")
        for i in self.finalSubstitutions:
            f.write(">S")
            f.write(str(i[0]))
            f.write(" ")
            f.write(str(i[1]))
            f.write(" ")
            f.write(str(i[2]))
            f.write('\n')
        for i in self.finalInsertions:
            f.write(">I")
            f.write(str(i[0]))
            f.write(" ")
            f.write(str(i[1]))
            f.write('\n')
        for i in self.finalDeletions:
            f.write(">D")
            f.write(str(i[0]))
            f.write(" ")
            f.write(str(i[1]))
            f.write('\n')



referenceGenome = "sample_1000\sample_1000_reference_genome.fasta"
# referenceGenome = "project1a_10000_reference_genome.fasta"

# reads = "sample_1000\sample_1000_no_error_single_reads.fasta"
# reads = "sample_1000\sample_1000_with_error_single_reads.fasta"
# reads = "sample_1000\sample_1000_with_error_paired_reads.fasta"
# reads = "sample_1000\sample_1000_no_error_paired_reads REDUCED.fasta"
# reads = "sample_1000\sample_1000_no_error_paired_reads.fasta"
# reads = "sample_1000\TEST.fasta"
# reads = "sample_1000\DELETION.fasta"
reads = "sample_1000\SINGLE_DELETION.fasta"
# reads = "sample_1000\INSERTION.fasta"
# reads = "sample_1000\REDUCED.fasta"
# reads = "project1a_10000_with_error_paired_reads.fasta"
# reads = "test/project1a_10000_with_error_paired_reads_TEST.fasta"

if len(sys.argv) > 1:
    referenceGenome = sys.argv[1]
    reads = sys.argv[2]

aligner = readAligner(genomeFile=referenceGenome, readFile=reads)