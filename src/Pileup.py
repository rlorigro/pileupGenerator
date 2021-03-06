import pysam
from pyfaidx import Fasta
from collections import defaultdict
from datetime import datetime
from PIL import Image
import sys
import numpy
import copy
"""
This  module takes an alignment file and produces a pileup across all alignments in a query region, and encodes the
pileup as an image where each x pixel corresponds to base position and y corresponds to coverage depth
"""

class Pileup:
    '''
    A child of PileupGenerator which contains all the information and methods to produce a pileup at a given
    position, using the pysam and pyfaidx object provided by PileupGenerator
    '''

    def __init__(self,sam,fasta,chromosome,queryStart,flankLength,outputFilename,label,insertLengths,deleteLengths,coverageCutoff,mapQualityCutoff,windowCutoff):
        self.length = flankLength*2+1
        self.label = label
        # self.outputLabel = label
        self.insertLengths = insertLengths
        self.deleteLengths = deleteLengths
        self.coverageCutoff = coverageCutoff
        self.outputFilename = outputFilename
        self.queryStart = queryStart
        self.queryEnd = queryStart + self.length -1
        self.chromosome = chromosome
        self.mapQualityCutoff = mapQualityCutoff
        self.windowCutoff = windowCutoff

        # pysam uses 0-based coordinates
        self.localReads = sam.fetch("chr"+self.chromosome, start=self.queryStart, end=self.queryEnd)

        # pyfaidx uses 1-based coordinates
        self.refSequence = fasta.get_seq(name="chr"+self.chromosome, start=self.queryStart, end=self.queryEnd)
        self.outputRefSequence = copy.deepcopy(self.refSequence)
        self.referenceRGB = list()

        # stored during cigar string parsing to save time
        self.inserts = defaultdict(list)

        self.deltaRef  = [1,0,1,0,0,0,0,0,0,0,0]    # key for whether reference advances
        self.deltaRead = [1,1,0,0,0,0,0,0,0,0,0]    # key for whether read sequence advances
                                                    #  ['M','I','D','N','S','H','P','=','X','B','NM']

        self.noneChar = '_'       # character to use for empty positions in the pileup
        self.noneLabel = '0'      # character to use for (non variant called) inserts in the label

        self.SNPtoRGB = {'M': [255,255,255],
                         'A': [255,0,  0],
                         'C': [255,255,0],
                         'G': [0,  255,0],
                         'T': [0,  0,  255],
                         'I': [255,0,  255],
                         'D': [0,  255,255],
                         'N': [0,  0,  0],  # redundant in case of read containing 'N'... should this be independent?
               self.noneChar: [0,  0,  0],}

        self.sortingKey = {'M':5,
                           'A':0,
                           'C':1,
                           'G':2,
                           'T':3,
                           'I':6,
                           'D':4,
                           'N':7}

        self.RGBtoSNP = [[['N', 'T'], ['G', 'D']], [['A', 'I'], ['C', 'M']]]

        self.RGBtoSNPtext = [[[' ', 'T'], ['G', 'D']], [['A', '_'], ['C', '.']]]

        self.pileupColumns = dict()
        self.insertColumns = dict()
        self.pileupImage = list()

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


    def addPileupEntry(self,i,pileupEncoding,quality):
        if i not in self.pileupColumns:
            self.initializePileupColumn(i)

        self.pileupColumns[i][pileupEncoding].append(tuple(self.SNPtoRGB[pileupEncoding]+[quality]))


    def addInsertEntry(self,i,readCharacters,n,qualities):
        '''
        Given a column index (i), the read identity at that position (read characters), their length (n) and qualities,
        update the column's pool of sorted segments (M,A,C,G,T,N)
        :param i:
        :param readCharacters:
        :param n:
        :param qualities:
        :return:
        '''
        if i not in self.insertColumns:
            self.insertColumns[i] = list()

        for k in range(n):
            if k >= len(self.insertColumns[i]):     # list is empty
                self.insertColumns[i].append({'A':list(),
                                              'C':list(),
                                              'G':list(),
                                              'T':list(),
                                              'N':list()})

            self.insertColumns[i][k][readCharacters[k]].append(tuple(self.SNPtoRGB[readCharacters[k]]+[qualities[k]]))


    def initializePileupColumn(self,i):
        self.pileupColumns[i] = {'M':list(),
                                 'D':list(),
                                 'A':list(),
                                 'C':list(),
                                 'G':list(),
                                 'T':list(),
                                 'N':list()}


    def getPileupEncoding(self,cigarCode, refCharacter, readCharacter):
        if cigarCode == 'D':
            return 'D'
        elif readCharacter == refCharacter:
            return 'M'
        else:
            return readCharacter


    def generateRGBtoSNP(self):
        '''
        Generate inverse of SNPtoRGB automatically, from any SNPtoRGB dictionary
        :return:
        '''

        self.RGBtoSNP = [[[None for i in range(2)] for j in range(2)] for k in range(2)]

        for item in self.SNPtoRGB.items():
            bits = [int(i/255) for i in item[1]]

            i1,i2,i3 = bits
            self.RGBtoSNP[i1][i2][i3] = item[0]


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
        quality = (1-(10 ** ((quality)/-10)))*255   # convert to probability and scale to 0-255 for pixel encoding
        quality = int(round(quality))

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

                            self.addInsertEntry(self.absolutePosition-1, readCharacters, n, qualities)

                        self.readPosition += self.deltaRead[snp]*n

                    elif snp < 4:   # not an insert
                        for i in range(n):

                            if self.absolutePosition >= self.queryStart and self.absolutePosition <= self.queryEnd:
                                readCharacter = self.readSequence[self.readPosition]
                                refCharacter = self.refSequence[self.relativeIndexRef]
                                cigarCode = self.cigarLegend[snp]
                                baseQuality = readQualities[self.readPosition]

                                quality = self.calculateQuality(baseQuality, mapQuality)

                                pileupEncoding = self.getPileupEncoding(cigarCode,refCharacter,readCharacter)

                                self.addPileupEntry(self.absolutePosition, pileupEncoding, quality)

                                self.relativeIndexRef += self.deltaRef[snp]

                                self.refPosition += self.deltaRef[snp]
                                self.readPosition += self.deltaRead[snp]
                            else:
                                self.refPosition += self.deltaRef[snp]
                                self.readPosition += self.deltaRead[snp]

                            self.absolutePosition = self.refPosition+self.refStart

                else:  # stop iterating after query region
                    break



    def encodeColumn(self,columnIndex,cigarSegments):
        '''
        Take a column-specific dictionary and turn it into a RGBA pixel vector
        '''

        if columnIndex is None:
            ntRef = 'I'
            fillerPixel = tuple(self.SNPtoRGB['I']+[255])
        else:
            ntRef = self.refSequence[columnIndex]
            fillerPixel = tuple(self.SNPtoRGB[self.noneChar]+[255])

        column = [tuple(self.SNPtoRGB[ntRef]+[255])]

        length = 0
        for item in sorted(cigarSegments.items(), key=lambda x: self.sortingKey[x[0]]):
            cigarCode = item[0]
            segment = item[1]

            column += segment
            length += len(segment)

        column += [fillerPixel]*(self.coverageCutoff-length)
        self.pileupImage.append(column)


    def generatePileupImage(self):
        '''
        Given that all reads have been iterated through and their cigar codes have been stored in segments, create a
        2D list with the width of the coverage cutoff parameter and length of the window cutoff parameter.
        Additionally, identify insert and delete positions (stored during iteration) within the pileup and update the
        image + label as needed.
        :return:
        '''

        offset = 0

        for entry in sorted(self.pileupColumns.items()):
            absolutePosition = entry[0]                             # chromosome position
            cigarSegments = entry[1]                                # dictionary of counts for M,D,A,C,G,T
            columnIndex = absolutePosition - self.queryStart        # position in the query window

            self.encodeColumn(columnIndex,cigarSegments)

            if columnIndex in self.deleteLengths:                   # update delete labels
                label = self.label[columnIndex+offset]
                length = self.deleteLengths[columnIndex]

                for c in range(1,length+1):
                    self.label = self.label[:columnIndex+c+offset] + label + self.label[columnIndex+c+offset+1:]

            if absolutePosition in self.insertColumns:              # update insert labels
                n = len(self.insertColumns[absolutePosition])
                self.labelInsertRegion(columnIndex,offset,n)

                for i,entry in enumerate(self.insertColumns[absolutePosition]):     # n insert columns
                    offset += 1
                    self.encodeColumn(None,entry)


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


    def savePileupRGB(self,filename):
        '''
        Save the pileup binary array as a bitmap image using a gray color map
        :param filename:
        :return:
        '''

        if not filename.endswith(".png"):
            filename += ".png"

        self.pileupImage = self.pileupImage[:self.windowCutoff]
        self.pileupImage = [row[:self.coverageCutoff] for row in self.pileupImage]

        image = Image.new("RGBA",(self.coverageCutoff,self.windowCutoff))
        pixels = image.load()

        jlength = len(self.pileupImage)

        for i in range(image.size[0]):
            for j in range(image.size[1]):

                if j < jlength:
                    pixels[i,j] = self.pileupImage[j][i] if i < len(self.pileupImage[j]) else tuple(self.SNPtoRGB[self.noneChar]+[255])
                else:
                    pixels[i,j] = tuple(self.SNPtoRGB[self.noneChar]+[255])
        image.save(filename,"PNG")


    def RGBtoBinary(self,rgb):
        return [int(value/255) for value in rgb]


    def RGBtoSortingKey(self,rgb):
        i1,i2,i3 = self.RGBtoBinary(rgb)
        code = self.RGBtoSNP[i1][i2][i3]
        return self.sortingKey[code]


    def getOutputLabel(self):
        blankLength = self.windowCutoff - len(self.label)

        self.label += self.noneLabel*blankLength
        return self.label


    def decodeRGB(self,filename):
        '''
        Read a RGB and convert to a text alignment
        :param filename:
        :return:
        '''

        img = Image.open(filename)          # <-- switch this to RGBA
        pixels = numpy.array(img.getdata())
        text = list()

        width,height = img.size
        depth = 4

        pixels = numpy.reshape(pixels,(height,width,depth))

        for h in range(height):
            row = ''
            for w in range(width):
                r,g,b,a = self.RGBtoBinary(pixels[h][w])

                row += self.RGBtoSNPtext[r][g][b]
            text.append(row)

        return text


class PileUpGenerator:
    '''
    Creates pileups of aligned reads given a SAM/BAM alignment and FASTA reference
    '''

    def __init__(self,alignmentFile,referenceFile):
        self.sam = pysam.AlignmentFile(alignmentFile,"rb")
        self.fasta = Fasta(referenceFile,as_raw=True,sequence_always_upper=True)


    def generatePileup(self,chromosome,position,flankLength,outputFilename,label,insertLengths,deleteLengths,coverageCutoff,mapQualityCutoff,windowCutoff):
        '''
        Generate a pileup at a given position
        :param queryStart:
        :param regionLength:
        :return:
        '''

        queryStart = position-flankLength +1

        chromosome = str(chromosome)


        # startTime = datetime.now()
        pileup = Pileup(self.sam,self.fasta,chromosome,queryStart,flankLength,outputFilename,label,insertLengths=insertLengths,deleteLengths=deleteLengths,windowCutoff=windowCutoff,coverageCutoff=coverageCutoff,mapQualityCutoff=mapQualityCutoff)


        # print(datetime.now() - startTime, "initialized")
        pileup.iterateReads()
        # print(datetime.now() - startTime, "drafted")
        pileup.generatePileupImage()
        # print(datetime.now() - startTime, "finalized")
        pileup.savePileupRGB(outputFilename)
        # print(datetime.now() - startTime, "encoded and saved")
        # print()

        # ----- UNCOMMENT FOR TEXT DECODING OF IMAGES ------
        # print(outputFilename)
        # label = pileup.getOutputLabel()
        #
        # rows = pileup.decodeRGB(outputFilename + ".png")
        # for r,row in enumerate(rows):
        #     print(label[r],row)
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


