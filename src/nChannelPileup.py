import pysam
from pyfaidx import Fasta
from collections import defaultdict
from datetime import datetime
import sys
import numpy
import copy
import csv
from scipy import sparse
from scipy import misc

"""
This  module takes an alignment file and produces a pileup across all alignments in a query region, and encodes the
pileup as an image where each x pixel corresponds to base position and y corresponds to coverage depth
"""

class Pileup:
    '''
    A child of PileupGenerator which contains all the information and methods to produce a pileup at a given
    position, using the pysam and pyfaidx object provided by PileupGenerator
    '''

    def __init__(self,sam,fasta,chromosome,queryStart,variantPosition,flankLength,outputFilename,coverageCutoff,mapQualityCutoff,windowCutoff,coverageThreshold):
        self.length = flankLength*2+1
        self.flankLength = flankLength
        self.coverageCutoff = coverageCutoff
        self.outputFilename = outputFilename
        self.variantPosition = variantPosition
        self.queryStart = queryStart
        self.queryEnd = queryStart + self.length - 1
        self.chromosome = chromosome
        self.mapQualityCutoff = mapQualityCutoff
        self.windowCutoff = windowCutoff
        self.coverageThreshold = coverageThreshold

        # pysam fetch reads
        self.localReads = sam.fetch("chr"+self.chromosome, start=self.queryStart, end=self.queryEnd)

        # pyfaidx fetch reference sequence for query region
        self.refSequence = fasta.get_seq(name="chr"+self.chromosome, start=self.queryStart, end=self.queryEnd)
        self.outputRefSequence = copy.deepcopy(self.refSequence)
        self.referenceRGB = list()

        self.deltaRef  = [1,0,1,0,0,0,0,0,0,0,0]    # key for whether reference advances
        self.deltaRead = [1,1,0,0,0,0,0,0,0,0,0]    # key for whether read sequence advances
                                                    #  ['M','I','D','N','S','H','P','=','X','B','NM']

        self.insertChar = '-'
        self.noneChar = '_'       # character to use for empty positions in the text pileup

        self.SNPtoRGB = {'M': [1.0,0.0,0.0,0.0,0.0,0.0,0.0],  # Match channel, 1=mismatch 0=match
                         'A': [0.0,1.0,0.0,0.0,0.0,0.0,0.0],
                         'C': [0.0,0.0,1.0,0.0,0.0,0.0,0.0],
                         'G': [0.0,0.0,0.0,1.0,0.0,0.0,0.0],
                         'T': [0.0,0.0,0.0,0.0,1.0,0.0,0.0],
                         'I': [0.0,0.0,0.0,0.0,0.0,1.0,0.0],
                         'D': [0.0,0.0,0.0,0.0,0.0,0.0,1.0],
                         'N': [0.0,0.0,0.0,0.0,0.0,0.0,1.0],  # redundant in case of read containing 'N'... should this be independent?
               self.noneChar: [0.0,0.0,0.0,0.0,0.0,0.0,0.0]}

        self.emptyChannelVector = [0.0 for i in range(len(self.SNPtoRGB['M']))]

        self.channelKey = {'M':0,
                           'A':1,
                           'C':2,
                           'G':3,
                           'T':4,
                           'I':5,
                           'D':6,
                           'N':6,
                 self.noneChar:None}

        self.decodeMap = ['M', 'A', 'C', 'G', 'T', self.insertChar, 'D', self.noneChar]
        self.decodeIndex = list(range(len(self.decodeMap)))

        self.noneAlpha = [1]
        self.insertAlpha = [1]
        self.referenceAlpha = [1]

        self.inserts = dict()
        self.pileupImage = [[self.SNPtoRGB[self.noneChar]+self.noneAlpha]*(self.windowCutoff) for j in range(self.coverageCutoff)]  # mixing tuples (below) and lists... bad!
        self.cigarLegend = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B', '?']  # ?=NM

        self.cigarTuples = None
        self.refStart = None
        self.refEnd = None
        self.refPositions = None
        self.refPosition = None
        self.readSequence = None
        self.readPosition = None
        self.relativeIndex = None
        self.relativeIndexRef = None
        self.pileupEnds = dict()
        self.packMap = dict()
        self.refAnchors = list()
        self.coveragePerColumn = defaultdict(int)
        self.arrayShape = None
        self.nInsertsBeforeVariant = 0


    def generateDecodeMap(self):
        '''
        Generate inverse of SNPtoRGB automatically, from any SNPtoRGB dictionary
        :return:
        '''

        decodeMap = [None for n in range(len(self.SNPtoRGB.keys()))]

        for key in self.SNPtoRGB:
            vector = self.SNPtoRGB[key]
            i = vector.index(1)

            decodeMap[i] = key


    def generateReadMapping(self,r,startIndex):
        '''
        Find the topmost available space to insert into the pileup
        :param r:
        :param startIndex:
        :return:
        '''

        i = 0
        unmapped = True
        for item in sorted(self.pileupEnds.items(), key=lambda x: x[1]):
            i += 1
            end = item[1]

            if end < startIndex:
                self.packMap[r] = item[0]   # store for use during insert homogenization
                unmapped = False
                break

        if unmapped:
            self.packMap[r] = i


    def addPileupEntry(self, r, snp, n, columnIndex, quality):
        '''
        For a given read and SNP, add the corresponding encoding to the pileup array
        :param r:
        :param snp:
        '''

        # if self.packMap[r] == 0:
        #     print(snp,self.relativeIndexRef,self.absolutePosition, self.variantPosition)

        if snp < 4 and snp != 1:
            index = self.relativeIndexRef
            encoding = None
            coverage = None
            self.pileupEnds[self.packMap[r]] = index

            if snp == 0:                                            # match
                nt = self.readSequence[self.readPosition]
                ntRef = self.refSequence[self.relativeIndexRef]

                # if self.packMap[r] == 0:
                #     print(ntRef)

                if nt != ntRef and nt != 'N':     #mismatch AND NOT A SEQUENCING ERROR
                    encoding = list(self.SNPtoRGB['M'])         # Mismatch (specify the alt)
                else:
                    encoding = list(self.emptyChannelVector)

                ntChannel = self.channelKey[nt]

                encoding[ntChannel] = 1.0

                coverage = 1

            elif snp == 2:                                          # delete
                encoding = list(self.SNPtoRGB['D'])
                coverage = 1

                encoding[self.channelKey['M']] = 1.0    # set mismatch channel to 1

            elif snp == 3:                                          # refskip
                encoding = list(self.SNPtoRGB['N'])
                coverage = 0

            self.coveragePerColumn[index] += coverage

            encoding = encoding+[quality]  # append the quality Alpha value and store as tuple

            self.pileupImage[self.packMap[r]][index] = encoding      # Finally add the code to the pileup


    def addInsertEntry(self, readCharacters,n,qualities, r):
        '''
        Given a column index (i), create a dictionary of (read index):characters
        :param i:
        :param readCharacters:
        :param n:
        :param qualities:
        :return:
        '''
        i = self.relativeIndexRef
        if i not in self.inserts:
            self.inserts[i] = []

        for c, character in enumerate(readCharacters):
            if c >= len(self.inserts[i]):   # no insert column exists yet
                self.inserts[i] = self.inserts[i] + [[self.SNPtoRGB['I']+self.insertAlpha]*(self.coverageCutoff+1)]

            encoding = list(self.SNPtoRGB[character])
            encoding[self.channelKey['M']] = 1.0        # set mismatch channel to 1

            channelIndex = self.channelKey['I']

            encoding[channelIndex] = 1.0
            self.inserts[i][c][self.packMap[r]+1] = encoding + [qualities[c]]


    def getPileupEncoding(self,cigarCode, refCharacter, readCharacter):
        if cigarCode == 'D':
            return 'D'
        elif readCharacter == refCharacter:
            return 'M'
        else:
            return readCharacter


    def iterateReads(self):
        '''
        For all the reads mapped to the region specified, (subsample and) collect their alignment data and build a draft
        pileup
        '''

        pileupIteratorIndex = 0

        for r, read in enumerate(self.localReads):
            if read.mapping_quality < self.mapQualityCutoff:
                continue
            self.parseRead(pileupIteratorIndex, read)
            pileupIteratorIndex += 1


    def calculateQuality(self,baseQuality,mapQuality):
        '''
        Approximate the overall quality of the given read character based on its mapping and sequence quality
        :param baseQuality:
        :param mapQuality:
        :return:
        '''
        quality = min(baseQuality,mapQuality)       # calculate minimum of quality for base and read quality
        quality = (1-(10 ** ((quality)/-10)))   # convert to probability and scale to 0-255 for pixel encoding

        return quality


    def parseRead(self,r,read):
        '''
        Create a draft pileup representation of a read using cigar string and read data. This method iterates through
        a read and its cigar string from the start (so it needs to be optimized for long reads).
        :param r:
        :param read:
        :return:
        '''

        refPositions = read.get_reference_positions()

        readQualities = read.query_qualities
        mapQuality = read.mapping_quality

        firstPass = True

        if len(refPositions) == 0:
            # sys.stderr.write("WARNING: read contains no reference alignments: %s\n" % read.query_name)
            pass
        else:
            self.refStart = refPositions[0] + 1
            self.cigarTuples = read.cigartuples

            refPositions = None  # dump

            self.readSequence = read.query_alignment_sequence

            self.refPosition = 0        # ???
            self.readPosition = 0       # index of the queried read sequence
            self.relativeIndexRef = 0   # index of the reference sequence to compare mismatches with

            if self.refStart >= self.queryStart:                         # if the read starts during the query region
                self.relativeIndexRef = self.refStart-self.queryStart    # add the difference to ref index

            for c, entry in enumerate(self.cigarTuples):
                snp = entry[0]
                n = entry[1]

                self.absolutePosition = self.refPosition+self.refStart
                if self.absolutePosition < self.queryStart-n:           # skip by blocks if not in query region yet
                    self.refPosition += self.deltaRef[snp]*n
                    self.readPosition += self.deltaRead[snp]*n

                elif self.absolutePosition >= self.queryStart-n and self.absolutePosition < self.queryEnd:    # read one position at a time for query region
                    if snp == 1:    # insert
                        if self.absolutePosition > self.queryStart:
                            readCharacters = self.readSequence[self.readPosition:self.readPosition+n]
                            baseQualities = readQualities[self.readPosition:self.readPosition+n]

                            qualities = [self.calculateQuality(baseQuality, mapQuality) for baseQuality in baseQualities]

                            self.addInsertEntry(readCharacters, n, qualities, r)

                        self.readPosition += self.deltaRead[snp]*n


                    elif snp < 4:   # not an insert
                        for i in range(n):

                            if self.absolutePosition >= self.queryStart and self.absolutePosition <= self.queryEnd:
                                if firstPass:
                                    firstPass = False
                                    self.generateReadMapping(r, self.relativeIndexRef)

                                    if self.packMap[r] >= self.coverageCutoff:   # skip unmapped read if exceeds coverage
                                        return

                                baseQuality = readQualities[self.readPosition]

                                quality = self.calculateQuality(baseQuality, mapQuality)

                                self.addPileupEntry(r, snp, n, self.relativeIndexRef, quality)

                                self.relativeIndexRef += self.deltaRef[snp]
                                self.refPosition += self.deltaRef[snp]
                                self.readPosition += self.deltaRead[snp]
                            else:
                                self.refPosition += self.deltaRef[snp]
                                self.readPosition += self.deltaRead[snp]

                            self.absolutePosition = self.refPosition+self.refStart

                else:  # stop iterating after query region
                    break


    def encodeReference(self):
        '''
        Encode the reference sequence into RGB triplets and add add it as the header line to the pileup RGB
        :return:
        '''

        for c,character in enumerate(self.refSequence):
            encoding = self.SNPtoRGB[character]+self.referenceAlpha
            self.referenceRGB.append(encoding)


    def cleanInsertColumns(self,i):
        '''
        If an insert column position doesn't have read data to its left, it should be a None,
        :param i:
        :return:
        '''

        for c in range(len(self.inserts[i])):
            for r in range(self.coverageCutoff+1):
                prevEntry = self.pileupImage[r][i-1][:-1]   # exclude quality

                try:
                    index = prevEntry.index(1)
                except(ValueError):
                    index = len(self.decodeMap) - 1

                prevEntry = self.decodeMap[index]

                if prevEntry == self.noneChar: # and insertEntry not in charSet:  # previous entry decodes to None character AND insert entry is not a true sequence
                    self.inserts[i][c][r] = self.SNPtoRGB[self.noneChar]+self.noneAlpha


    def savePileupRGB(self,filename):
        '''
        Save the pileup binary array as a bitmap image using a gray color map
        :param filename:
        :return:
        '''

        if not filename.endswith(".png"):
            filename += ".png"

        self.pileupImage = self.pileupImage[:self.coverageCutoff]
        self.pileupImage = [row[:self.windowCutoff] for row in self.pileupImage]

        self.encodeReference()
        self.pileupImage = [self.referenceRGB] + self.pileupImage

        pileupArray = [[None for m in range(self.windowCutoff)] for n in range(self.coverageCutoff)]

        # count insert columns that appear before variant column
        for i in range(0,self.flankLength):
            if i in self.inserts:
                self.nInsertsBeforeVariant += len(self.inserts[i])

        # print(self.nInsertsBeforeVariant)

        image_iterator = -self.nInsertsBeforeVariant
        i = 0
        while image_iterator < self.windowCutoff:
            if i in self.inserts:
                insert_rows = self.inserts[i]
                self.cleanInsertColumns(i)

                for insert in insert_rows:
                    for j in range(self.coverageCutoff):
                        if image_iterator >= 0:
                            pileupArray[j][image_iterator] = insert[j] if j < len(insert) else self.SNPtoRGB[self.noneChar]+self.noneAlpha

                    image_iterator += 1
                    self.refAnchors.append(i-1)
                    if image_iterator >= self.windowCutoff:
                        break

            if image_iterator >= self.windowCutoff:
                break

            for j in range(self.coverageCutoff):
                if image_iterator >= 0:
                    if j < self.coverageCutoff:
                        pileupArray[j][image_iterator] = self.pileupImage[j][i] if i < len(self.pileupImage[j]) else self.SNPtoRGB[self.noneChar]+self.noneAlpha
                    else:
                        pileupArray[j][image_iterator] = self.SNPtoRGB[self.noneChar]+self.noneAlpha

            if i < len(self.refSequence):
                self.refAnchors.append(i)
            else:
                self.refAnchors.append(len(self.refSequence)-1)

            image_iterator += 1
            i += 1


        pileupArray = numpy.array(pileupArray)

        # print(pileupArray.shape)

        #------
        pileupArray2d = pileupArray.reshape((pileupArray.shape[0],-1))
        self.arrayShape = pileupArray.shape

        misc.imsave(self.outputFilename+".png",pileupArray2d,format="PNG")
        #------

        #------
        # pileupArray2d = pileupArray.reshape((pileupArray.shape[0],-1))
        # sparseMatrix = sparse.csc_matrix(pileupArray2d)
        #
        # arrayFilename = self.outputFilename + ".npz"
        # sparse.save_npz(arrayFilename,sparseMatrix)
        #------

        # numpy.save(self.outputFilename,pileupArray) # <- waste of space


class PileUpGenerator:
    '''
    Creates pileups of aligned reads given a SAM/BAM alignment and FASTA reference
    '''

    def __init__(self,alignmentFile, referenceFile):    #, smry_ref_pos_file_writer):
        self.sam = pysam.AlignmentFile(alignmentFile,"rb")
        self.fasta = Fasta(referenceFile,as_raw=True,sequence_always_upper=True)

    def generatePileup(self,chromosome,
                       position,
                       flankLength,
                       outputFilename,
                       coverageCutoff,
                       mapQualityCutoff,
                       windowCutoff,
                       coverageThreshold):
        '''
        Generate a pileup at a given position
        :param queryStart:
        :param regionLength:
        :return:
        '''

        queryStart = position-flankLength + 1

        chromosome = str(chromosome)


        startTime = datetime.now()
        pileup = Pileup(self.sam,
                        self.fasta,
                        chromosome=chromosome,
                        queryStart=queryStart,
                        variantPosition=position+1,
                        flankLength=flankLength,
                        outputFilename=outputFilename,
                        windowCutoff=windowCutoff,
                        coverageCutoff=coverageCutoff,
                        mapQualityCutoff=mapQualityCutoff,
                        coverageThreshold=coverageThreshold)


        # print(datetime.now() - startTime, "initialized")
        pileup.iterateReads()
        # pileup.generateDecodeMap()
        # print(datetime.now() - startTime, "finalized")
        pileup.savePileupRGB(outputFilename)
        # print(datetime.now() - startTime, "encoded and saved")
        # print()

        return pileup.queryStart, pileup.refAnchors, pileup.arrayShape


