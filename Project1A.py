import sys
from copy import deepcopy
import copy
import time

class readAligner():
    kmerSize = 10
    genomeFile = None
    readFile = None
    reads = None
    genome = ""
    genomeMap = None
    allowedErrors = 3
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
    stepSize = 10
    currentRead = None

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
        self.genomeMap = genomeMap
        print("genome map complete")

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
        # print(genomeString)
        # print(readString)
        if len(genomeString) != len(readString):
            # print(genomeString)
            # print(readString)
            raise ValueError("Input strings must have the same length")
        mutation = []
        for i in range(len(genomeString)):
            if genomeString[i] != readString[i] and genomeString[i] != ' ':
                if (i+alignment, readString[i], genomeString[i]) not in self.substitutions:
                    # if i+alignment == 333012:
                    #     # print(self.currentRead)
                    self.substitutions[(i+alignment, readString[i], genomeString[i])] = 1
                else:
                    self.substitutions[(i+alignment, readString[i], genomeString[i])] +=1
                mutation.append((i+alignment, readString[i], genomeString[i]))
        return mutation

    def getStreak(self, positions):
        avg = 0
        avgCounter = 0
        # print(positions, "pos in streak")
        # print(positions == [[-10]]))
        # for i in positions:
        #     if i != [-10] and len(i) == 1:
        #         avg += i[0]
        #         avgCounter+=1

        # print(positions, "positions")
        positionsCopy = deepcopy(positions)
        # if avgCounter == 0:
            # print("in")
        for count2, subArray in enumerate(positions):
            if subArray == [-10]:
                continue
            for subArrayItem in subArray:
                counter = 0
                for count, subArray2 in enumerate(positions):
                    for j in range(-1, 2):
                        # print(subArrayItem+(count-count2)*self.stepSize + j, subArray2)
                        if (subArrayItem+(count-count2)*self.stepSize + j in subArray2):
                            # print("add")
                            counter+=1
                # print(subArrayItem, counter)
                if counter == 1 and len(subArray) > 1:
                    positionsCopy[count2].remove(subArrayItem)
            if len(positionsCopy[count2]) == 0:
                positionsCopy[count2] = [-10]
            pass
        # print(positionsCopy)
        positions = deepcopy(positionsCopy)
        for i in positions:
            if i != [-10] and len(i) == 1:
                avg += i[0]
                avgCounter+=1



        if positions == [[-10]] or avgCounter == 0:
            # print(positions)
            # print("skippppp")
            return []
        avg /= avgCounter
        # print(avg, "avg")
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

        # print(positionsReduced, "reduced")

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
                if streak == False: # start new streak
                    # print(i)
                    streak = True
                    currentStreakStartIndex = count
                    currentStreakSum = i
                    currentStreakCount = 1

                    streakIndexes.append(currentStreakStartIndex)
                    streakCounts.append(currentStreakCount)
                    streakAvgs.append(currentStreakSum/currentStreakCount)

                elif streak == True:
                    if abs(i - last)/self.stepSize < 2: #account for mutation
                        currentStreakSum += i
                        currentStreakCount += 1
                        if count == (len(positionsReduced)-1):
                            streak = False
                            streakAvgs[-1] = (currentStreakSum/currentStreakCount)
                            streakCounts[-1] += 1
                    else: # current streak ended by new number
                        currentStreakStartIndex = count
                        currentStreakSum = i
                        currentStreakCount = 1


                        streakAvgs.append(currentStreakSum/currentStreakCount)
                        streakIndexes.append(currentStreakStartIndex)
                        streakCounts.append(currentStreakCount)
                last = i
            elif i == -10: #streak ended by -10
                if streak == True:
                    streak = False
                    # streakAvgs.append(currentStreakSum/currentStreakCount)
                    # streakIndexes.append(currentStreakStartIndex)
                    # streakCounts.append(currentStreakCount)


        # print(streakCounts, "streakCounts")
        # print(streakIndexes, "streakIndexes")

        maxStreakCountIndex = streakCounts.index(max(streakCounts))
        maxStreakStartIndex = streakIndexes[maxStreakCountIndex]

        # print(streakCounts)
        # print(streakIndexes)
        # print(maxStreakStartIndex)

        streakIndexes.pop(maxStreakCountIndex)
        streakCounts.pop(maxStreakCountIndex)
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
        readCount = 0
        keys = list(self.genomeMap.keys())[0:] #THIS LINE TAKES FOREVER
        for count, read in enumerate(self.reads):
            # print(lastReadPos)
            self.currentRead = read
            readSplit = read[0].split("/")
            # lastReadPos = 0
            if len(readSplit) > 1:
                if read[0].split("/")[1] == '1':
                    lastReadPos = 0
            else:
                lastReadPos = 0
            # lastReadPos = 771746

            currentTime = time.perf_counter() - startTime
            # print(read[0], '\n', "readPos:", lastReadPos, '\n', "currentTime:", round(currentTime, 3), '\n', "readTime:",  round(currentTime - lastCurrentTime, 3))
            readCount += 1

            if readCount % 10000 == 0:
                if count == len(self.reads) - 1:
                    print(read[0], round(currentTime, 3), round(currentTime - lastCurrentTime, 10), end='\n')
                else:
                    print(read[0], round(currentTime, 3), round(currentTime - lastCurrentTime, 10), end='\n')
                lastCurrentTime = currentTime


            insertionState = False
            deletionState = False

            positions = [[-10]] * (int((len(read[1])/1-self.kmerSize)/self.stepSize) + 1)
            for kmerPos in range(0, len(read[1])-self.kmerSize+1, self.stepSize):
                # keys = list(self.genomeMap.keys())[lastReadPos:]
                # keys = list(self.genomeMap.keys())[0:] #THIS LINE TAKES FOREVER
                counter = 0

                # for key in keys:
                #     if read[1][kmerPos:kmerPos+self.kmerSize] == key:
                #         currentTime = time.perf_counter() - startTime
                #         # print(read[1])
                #         lastCurrentTime = currentTime
                #         if positions[int(kmerPos/self.stepSize)] == [-10]:
                #             positions[int(kmerPos/self.stepSize)] = deepcopy(self.genomeMap[key])
                #         else:
                #             positions[int(kmerPos/self.stepSize)] += deepcopy(self.genomeMap[key])
                #         break

                kmerToCheck = read[1][kmerPos:kmerPos+self.kmerSize]

                if kmerToCheck in self.genomeMap:
                    if positions[int(kmerPos/self.stepSize)] == [-10]:
                        positions[int(kmerPos/self.stepSize)] = copy.copy(self.genomeMap[kmerToCheck])
                    else:
                        positions[int(kmerPos/self.stepSize)] += copy.copy(self.genomeMap[kmerToCheck])

            if len(positions) == 0:
                continue

            positionsFiltered = self.getStreak(positions)

            if positionsFiltered == []:
                continue

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
                    if abs(posItem-last) <= self.stepSize+1:
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
            
            for count, k in enumerate(finalPositions):
                if first == None:
                    if k != -10:
                        first = k
                    continue
                if k != first + 1*self.stepSize:
                    diff = (k - first)
                    if diff == 1+self.stepSize:
                        lastHamming = 999
                        bestDeletionLoc = None
                        bestDeletionIndex = 0
                        for i in range(k-self.stepSize, k+1):
                            # print(bestDeletionIndex)
                            deletionLoc = i
                            deletionIndex = int((deletionLoc-finalPositions[0]))

                            currentModifiedRead = read[1][:deletionIndex] + self.genome[deletionLoc] + read[1][deletionIndex:]
                            currentHamming = self.hammingDistance(currentModifiedRead, self.genome[finalPositionsTwo[0]:finalPositionsTwo[0]+len(currentModifiedRead)])
                            if currentHamming < lastHamming and currentHamming < 3:
                                lastHamming = currentHamming
                                bestDeletionLoc = deletionLoc
                                bestDeletionIndex = deletionIndex
                                modifiedRead = currentModifiedRead
                                if currentHamming == 0:
                                    break
                        if deletionLoc < 0:
                            print(read)
                            print(positions, "pos")
                            print(finalPositions, "final")
                            print(finalPositionsTwo, "finalTwo")
                            print("---")
                            break

                        if bestDeletionLoc != None:
                            deletion = (bestDeletionLoc, self.genome[finalPositions[0] + bestDeletionIndex])
                            if deletion in self.deletions:
                                self.deletions[deletion] += 1
                            else:
                                self.deletions[deletion] = 1
                            deletionState = True
                        else:
                            pass

                    elif diff == self.stepSize-1:
                        insertionLoc = k-1
                        lastHamming = 999
                        bestInsertionLoc = None
                        bestInsertionIndex = None
                        for i in range(k-self.stepSize, k+1):
                            insertionLoc = i
                            insertionIndex = int((insertionLoc-finalPositions[0]))+1
                            currentModifiedRead = read[1][:insertionIndex] + read[1][insertionIndex+1:]
                            currentHamming = self.hammingDistance(currentModifiedRead, self.genome[finalPositionsTwo[0]:finalPositionsTwo[0]+len(currentModifiedRead)])
                            if currentHamming == 0 or (currentHamming < lastHamming and currentHamming < 2):
                                lastHamming = currentHamming
                                bestInsertionLoc = insertionLoc
                                bestInsertionIndex = insertionIndex
                                modifiedRead = currentModifiedRead

                        if bestInsertionLoc != None:
                            insertion = (bestInsertionLoc, read[1][bestInsertionIndex])
                            # print(insertion)
                            if insertion in self.insertions:
                                self.insertions[insertion] += 1
                            else:
                                self.insertions[insertion] = 1
                            first = k
                            insertionState = True
                first = k

            if len(finalPositionsTwo) > 0:
                for pos in range(finalPositionsTwo[0], finalPositionsTwo[-1]+self.kmerSize):
                    # print("pos", pos)
                    self.coverageMap[pos] += 1

                lastReadPos = finalPositionsTwo[0]
                stringAlignment = 0
                genomeString = ""

                genomeString = self.genome[finalPositionsTwo[0]:finalPositionsTwo[0]+len(modifiedRead)]
                stringAlignment = self.hammingDistance(modifiedRead, genomeString)

                # stringAlignment = self.hammingDistance(modifiedRead, genomeString)
                if [-10] != finalPositionsTwo and stringAlignment <= self.allowedErrors:
                    # print("mutation")
                    self.getSubstitutions(modifiedRead, genomeString, finalPositionsTwo[0])
        
        # print("Deletions:" )
        # for i in self.deletions:
        #     print(i)
        #     print(i, self.deletions[i], self.coverageMap[i[0]])
        # # print("---")
        # # print("deletions: ", self.deletions)
        # print("Insertions: ")
        # for i in self.insertions:
        #     print(i, self.insertions[i], self.coverageMap[i[0]])
        # # print("---")
        # # print("Substitutions:")
        # print("Subsitutions")
        # for i in self.substitutions:
        #     print(i, self.substitutions[i], self.coverageMap[i[0]])
        
        for i in self.deletions:
            if self.coverageMap[i[0]] > 1 and self.deletions[i] > 1:
                self.finalDeletions[i] = self.deletions[i]
        for i in self.insertions:
            if self.coverageMap[i[0]] > 1 and self.insertions[i] > 1:
                self.finalInsertions[i] = self.insertions[i]
        for i in self.substitutions:
            if self.coverageMap[i[0]] == 0:
                self.coverageMap[i[0]] +=1


            if self.substitutions[i] / self.coverageMap[i[0]] >= 0.5 and self.coverageMap[i[0]] > 1 and self.substitutions[i] > 1:
                self.finalSubstitutions[i] = self.substitutions[i]


    def writeResults(self):
        f = open("mutationsRAW.txt", "w")
        for i in self.substitutions:
            f.write(">S")
            f.write(str(i[0]))
            f.write(" ")
            f.write(str(i[1]))
            f.write(" ")
            f.write(str(i[2]))
            f.write('\n')
        for i in self.insertions:
            f.write(">I")
            f.write(str(i[0]))
            f.write(" ")
            f.write(str(i[1]))
            f.write('\n')
        for i in self.deletions:
            f.write(">D")
            f.write(str(i[0]))
            f.write(" ")
            f.write(str(i[1]))
            f.write('\n')

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



# referenceGenome = "sample_1000\sample_1000_reference_genome.fasta"
# referenceGenome = "project1a_10000_reference_genome.fasta"
referenceGenome = "project1b_1000000_reference_genome.fasta"

# reads = "sample_1000\sample_1000_no_error_single_reads.fasta"
# reads = "sample_1000\sample_1000_with_error_single_reads.fasta"
# reads = "sample_1000\sample_1000_with_error_paired_reads.fasta"
# reads = "sample_1000\sample_1000_no_error_paired_reads REDUCED.fasta"
# reads = "sample_1000\sample_1000_no_error_paired_reads.fasta"
# reads = "sample_1000\TEST.fasta"
# reads = "sample_1000\DELETION.fasta"
# reads = "sample_1000\SINGLE_DELETION.fasta"
# reads = "sample_1000\INSERTION.fasta"
# reads = "sample_1000\REDUCED.fasta"
# reads = "project1a_10000_with_error_paired_reads.fasta"
# reads = "test/project1a_10000_with_error_paired_reads_TEST.fasta"
reads = "project1b_1000000_with_error_paired_reads.fasta"
# reads = "project1b_1000000_with_error_paired_reads_TEST.fasta"

if len(sys.argv) > 1:
    referenceGenome = sys.argv[1]
    reads = sys.argv[2]

aligner = readAligner(genomeFile=referenceGenome, readFile=reads)