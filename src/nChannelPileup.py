import pysam
from pyfaidx import Fasta
from collections import defaultdict
from datetime import datetime
from PIL import Image
import sys
import numpy
import copy
import csv
"""
This  module takes an alignment file and produces a pileup across all alignments in a query region, and encodes the
pileup as an image where each x pixel corresponds to base position and y corresponds to coverage depth
"""

class Pileup:
    '''
    A child of PileupGenerator which contains all the information and methods to produce a pileup at a given
    position, using the pysam and pyfaidx object provided by PileupGenerator
    '''

    def __init__(self,sam,fasta,smry_ref_pos_file_writer,chromosome,queryStart,flankLength,outputFilename,label,insertLengths,deleteLengths,coverageCutoff,mapQualityCutoff,windowCutoff):
        self.length = flankLength*2+1
        self.label = label
        self.insertLengths = insertLengths
        self.deleteLengths = deleteLengths
        self.coverageCutoff = coverageCutoff
        self.outputFilename = outputFilename
        self.queryStart = queryStart
        self.queryEnd = queryStart + self.length - 1
        self.chromosome = chromosome
        self.mapQualityCutoff = mapQualityCutoff
        self.windowCutoff = windowCutoff
        self.smry_ref_pos_file_writer = smry_ref_pos_file_writer

        # pysam fetch reads
        self.localReads = sam.fetch("chr"+self.chromosome, start=self.queryStart, end=self.queryEnd)

        # pyfaidx fetch reference sequence for query region
        self.refSequence = fasta.get_seq(name="chr"+self.chromosome, start=self.queryStart, end=self.queryEnd)
        self.outputRefSequence = copy.deepcopy(self.refSequence)
        self.referenceRGB = list()

        # stored during cigar string parsing to save time
        # self.inserts = defaultdict(list)

        self.deltaRef  = [1,0,1,0,0,0,0,0,0,0,0]    # key for whether reference advances
        self.deltaRead = [1,1,0,0,0,0,0,0,0,0,0]    # key for whether read sequence advances
                                                    #  ['M','I','D','N','S','H','P','=','X','B','NM']

        self.noneChar = '_'       # character to use for empty positions in the text pileup
        self.noneLabel = '0'      # character to use for (non variant called) inserts in the label

        self.SNPtoRGB = {'M': [1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
                         'A': [0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0],
                         'C': [0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0],
                         'G': [0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0],
                         'T': [0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0],
                         'I': [0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0],
                         'D': [0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0],
                         'N': [0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0],  # redundant in case of read containing 'N'... should this be independent?
               self.noneChar: [0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0],}

        self.decodeMap = ['M', 'A', 'C', 'G', 'T', '-', 'D', '_']
        self.decodeIndex = list(range(len(self.decodeMap)))

        self.noneAlpha = [1]
        self.insertAlpha = [1]
        self.referenceAlpha = [1]

        self.sortingKey = {'M':5,
                           'A':0,
                           'C':1,
                           'G':2,
                           'T':3,
                           'I':6,
                           'D':4,
                           'N':7}

        # self.pileupColumns = dict()
        self.inserts = dict()
        self.pileupImage = [[self.SNPtoRGB[self.noneChar]+self.noneAlpha]*(self.windowCutoff) for j in range(self.coverageCutoff)]  # mixing tuples (below) and lists... bad!
        # print(self.pileupImage)
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

    def generateDecodeMap(self):
        '''
        Generate inverse of SNPtoRGB automatically, from any SNPtoRGB dictionary
        :return:
        '''

        indexMap = list(range(len(self.SNPtoRGB.keys())))
        decodeMap = [None for n in range(len(self.SNPtoRGB.keys()))]
        print(indexMap)

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

        if snp < 4 and snp != 1:
            index = self.relativeIndexRef
            encoding = None
            self.pileupEnds[self.packMap[r]] = columnIndex

            if snp == 0:                                            # match
                nt = self.readSequence[self.readPosition]
                ntRef = self.refSequence[self.relativeIndexRef]

                if nt != ntRef:
                    encoding = self.SNPtoRGB[nt]  # Mismatch (specify the alt)
                else:
                    encoding = self.SNPtoRGB['M']


            elif snp == 2:                                          # delete
                encoding = self.SNPtoRGB['D']

            elif snp == 3:                                          # refskip
                encoding = self.SNPtoRGB['N']


            encoding = encoding+[quality]  # append the quality Alpha value and store as tuple

            # print('r',r,'packmap',self.packMap[r],'index',index,'cov cutoff',self.coverageCutoff,len(self.pileupImage),len(self.pileupImage[0]))
            # print(encoding)
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
        # print(self.relativeIndexRef, readCharacters, n, qualities, r)
        i = self.relativeIndexRef
        if i not in self.inserts:
            self.inserts[i] = []

        for c, character in enumerate(readCharacters):
            if c >= len(self.inserts[i]):   # no insert column exists yet
                self.inserts[i] = self.inserts[i] + [[self.SNPtoRGB['I']+self.insertAlpha]*(self.coverageCutoff+1)]


            # print("inserts",i,len(self.inserts))
            # print("insert column",c,len(self.inserts[i]))
            # print("qual",c,len(qualities))
            # print("packmap",self.packMap[r]+1,len(self.inserts[i][c]))
            # print("r",r)
            self.inserts[i][c][self.packMap[r]+1] = self.SNPtoRGB[character] + [qualities[c]]


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
        # quality = int(round(quality))

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


    def labelInsertRegion(self,c,offset,n):
        '''
        Label inserts as some combination of None labels and het/hom ref/alt depending on the length of the insert
        specified by the vcf.
        :param c:
        :param offset:
        :param n:
        :return:
        '''
        if c in self.insertLengths:                         # if the position is a variant site
            l = self.insertLengths[c]  # length of insert
            if l >= n:                                      # correctly modify the label to fit the insert
                labelInsert = self.label[c+offset]*n        # using the length of the called variant at pos.
            else:
                labelInsert = self.label[c+offset]*l+self.noneLabel*(n-l)

        else:
            labelInsert = self.noneLabel*n                  # otherwise the insert is labeled with None label

        self.label = self.label[:c+offset+1]+labelInsert+self.label[c+offset+1:]


    def encodeReference(self):
        '''
        Encode the reference sequence into RGB triplets and add add it as the header line to the pileup RGB
        :return:
        '''

        for character in self.refSequence:
            encoding = self.SNPtoRGB[character]+self.referenceAlpha
            self.referenceRGB.append(encoding)


    def cleanInsertColumns(self,i):
        print(self.inserts.keys())
        for c in range(len(self.inserts[i])):
            for r in range(self.coverageCutoff+1):
                prevEntry = self.pileupImage[r][i-1][:-1]   # exclude quality
                prevEntry = self.decodeMap[prevEntry.index(1)]
                insertEntry = self.inserts[i][c][r][:-1]    # exclude quality
                insertEntry = self.decodeMap[insertEntry.index(1)]

                print(i,r,prevEntry,insertEntry)
                if prevEntry == self.noneChar and insertEntry != ['I']:  # previous entry decodes to None character
                    self.inserts[i][c][r] = self.SNPtoRGB[self.noneChar]+self.noneAlpha

            # self.inserts[i].append(insertTemplate)
            # print(self.inserts[i])


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

        image_iterator = 0
        i = 0
        while image_iterator < self.windowCutoff:
            if i in self.deleteLengths:                     # update delete labels
                label = self.label[image_iterator]
                length = self.deleteLengths[i]

                for c in range(1, length+1):
                    self.label = self.label[:image_iterator+c] + label + self.label[image_iterator+c+1:]

            if i in self.inserts:
                insert_rows = self.inserts[i]
                n = len(self.inserts[i])
                self.labelInsertRegion(i-1,image_iterator-i,n)  # update insert labels

                self.cleanInsertColumns(i)

                for insert in insert_rows:
                    for j in range(self.coverageCutoff):
                        pileupArray[j][image_iterator] = insert[j] if j < len(insert) else self.SNPtoRGB[self.noneChar]+self.noneAlpha

                    image_iterator += 1
                    self.refAnchors.append(i)
                    if image_iterator >= self.windowCutoff:
                        break

            if image_iterator >= self.windowCutoff:
                break

            for j in range(self.coverageCutoff):
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

        numpy.save(self.outputFilename,pileupArray)
        imageFilename = self.outputFilename + ".npy"

        row = [imageFilename, self.queryStart] + self.refAnchors
        self.smry_ref_pos_file_writer.writerow(row)


    def RGBtoBinary(self,rgb):
        return [int(value/255) for value in rgb]

    def RGBtoSortingKey(self,rgb):
        i1,i2,i3 = self.RGBtoBinary(rgb)
        code = self.RGBtoSNP[i1][i2][i3]
        return self.sortingKey[code]

    def getOutputLabel(self):
        blankLength = self.windowCutoff - len(self.label) if len(self.label) < self.windowCutoff else 0

        # print(len(self.label),self.windowCutoff)
        label = self.label + self.noneLabel*blankLength
        label = label[:self.windowCutoff]
        # print(label)
        return label


    def decodeRGB(self,filename):
        '''
        Read a RGB and convert to a text alignment
        :param filename:
        :return:
        '''

        array = numpy.load(filename)

        width,height,depth = array.shape

        text = list()

        for w in range(width):
            row = list()
            for h in range(height):
                channels = array[w,h,:-1]
                index = int(numpy.sum(channels*self.decodeIndex))

                row.append(self.decodeMap[index])

            text.append(''.join(row))

            return text


class PileUpGenerator:
    '''
    Creates pileups of aligned reads given a SAM/BAM alignment and FASTA reference
    '''

    def __init__(self,alignmentFile, referenceFile, smry_ref_pos_file_writer):
        self.sam = pysam.AlignmentFile(alignmentFile,"rb")
        self.fasta = Fasta(referenceFile,as_raw=True,sequence_always_upper=True)
        self.smry_ref_pos_file_writer = smry_ref_pos_file_writer


    def generatePileup(self,chromosome,position,flankLength,outputFilename,label,insertLengths,deleteLengths,coverageCutoff,mapQualityCutoff,windowCutoff):
        '''
        Generate a pileup at a given position
        :param queryStart:
        :param regionLength:
        :return:
        '''

        queryStart = position-flankLength + 1

        chromosome = str(chromosome)


        startTime = datetime.now()
        pileup = Pileup(self.sam,self.fasta,self.smry_ref_pos_file_writer,chromosome,queryStart,flankLength,outputFilename,label,insertLengths=insertLengths,deleteLengths=deleteLengths,windowCutoff=windowCutoff,coverageCutoff=coverageCutoff,mapQualityCutoff=mapQualityCutoff)


        # print(datetime.now() - startTime, "initialized")
        pileup.iterateReads()
        # pileup.generateDecodeMap()
        # print(datetime.now() - startTime, "finalized")
        pileup.savePileupRGB(outputFilename)
        # print(datetime.now() - startTime, "encoded and saved")
        # print()

        # ----- UNCOMMENT FOR TEXT DECODING OF IMAGES ------
        # print(outputFilename)
        # label = pileup.getOutputLabel()
        #
        # rows = pileup.decodeRGB(outputFilename + ".npy")
        # print(label)
        # for r,row in enumerate(rows):
        #     print(row)
        # --------------------------------------------------

        return pileup.getOutputLabel()

#
# bamFile = "deePore/data/chr3_200k.bam"
# fastaFile = "deePore/data/chr3.fa"
#
# vcf_region = "3"
# window_size = 100
# filename = "deePore/data/test_packed/"
# labelString = '0'*(window_size+1)
# variantLengths = dict()
# positions =  [193657] # [77131,77174,66164,70895,77037,70217]


# for position in positions:
#     piler = PileUpGenerator(bamFile,fastaFile)
#     outputLabelString = piler.generatePileup(chromosome=vcf_region, position=position, flankLength=window_size,
#                                      outputFilename=filename, label=labelString, variantLengths=variantLengths,forceCoverage=True,coverageCutoff=100)

