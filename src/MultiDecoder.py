import numpy
from PIL import Image
import os

class Decoder:
    def __init__(self,summaryFilename,localPath=None):
        self.summaryFile = summaryFilename
        self.localPath = localPath
        if self.localPath[-1] != '/':
            self.localPath += '/'

        # matches
        self.decodeToSNPMap = ['M', 'A', 'C', 'G', 'T', 'I', 'D', '_']
        self.decodeToTextMap = ['.', 'A', 'C', 'G', 'T', '-', 'D', '_']
        self.decodeIndex = list(range(len(self.decodeToTextMap)))

        # no matches
        # self.decodeToSNPMap = ['A', 'C', 'G', 'T', 'I', 'D', '_']
        # self.decodeToTextMap = ['A', 'C', 'G', 'T', '-', 'D', '_']
        # self.decodeIndex = list(range(len(self.decodeToTextMap)))


        self.RGBtoText = [[[' ', 'T'], ['G', 'D']], [['A', '_'], ['C', '.']]]

        self.sortingKey = {'M':5,
                           'A':0,
                           'C':1,
                           'G':2,
                           'T':3,
                           'I':6,
                           'D':4,
                           'N':7}

        self.noneChar = '_'

        self.SNPtoRGB = {'M': [255,255,255],
                         'A': [255,0,  0],
                         'C': [255,255,0],
                         'G': [0,  255,0],
                         'T': [0,  0,  255],
                         'I': [127,127,127],
                         'D': [0,  255,255],
                         'N': [0,  255,255],  # redundant in case of read containing 'N'... should this be independent?
               self.noneChar: [0  ,0  ,0]}


    def decodeArraysToStdout(self):
        with open(self.summaryFile,'r') as file:
            for line in file:
                arrayPath, label = line.split(',')

                if os.path.exists(arrayPath):
                    pileupText = self.decodeArray(arrayPath)

                    print(arrayPath.split('/')[-1])
                    print(label)
                    for r, row in enumerate(pileupText):
                        print(row)
                else:
                    print("WARNING file not found: ",arrayPath)
                    pass


    def convertArraysToPNG(self):
        with open(self.summaryFile,'r') as file:
            for line in file:
                arrayPath, label = line.split(',')

                if os.path.exists(arrayPath):
                    self.convertArrayToRGBA(arrayPath)
                else:
                    print("WARNING file not found: ",arrayPath)
                    pass


    def decodeRGBstoStdout(self):
        with open(self.summaryFile,'r') as file:
            for line in file:
                imagePath, label = line.split(',')

                if os.path.exists(imagePath):
                    pileupText = self.decodeRGB(imagePath)

                    print(imagePath.split('/')[-1])
                    for r, row in enumerate(pileupText):
                        print(label[r], row)
                else:
                    print("WARNING file not found: ",imagePath)
                    pass


    def RGBtoBinary(self,rgb):
        return [int(value/255) for value in rgb]


    def RGBtoSortingKey(self,rgb):
        i1,i2,i3 = self.RGBtoBinary(rgb)
        code = self.RGBtoText[i1][i2][i3]
        return self.sortingKey[code]


    def convertArrayToRGBA(self,npyFilePath):
        array = numpy.load(npyFilePath)

        width,height,depth = array.shape

        image = Image.new("RGBA", (height, width))
        pixels = image.load()


        for w in range(width):
            for h in range(height):
                channels = array[w,h,:-1]
                quality = array[w,h,depth-1]     # q is always last channel
                decodeIndex = int(numpy.sum(channels*self.decodeIndex))
                snp = self.decodeToSNPMap[decodeIndex]

                rgb = self.SNPtoRGB[snp]
                rgba = rgb+[quality*255]
                # print(quality)
                rgba = tuple(map(int,rgba))
                pixels[h,w] = rgba

        pngFilePath = self.localPath + npyFilePath.split('/')[-1].split('.')[0] + ".png"
        print(pngFilePath)
        image.save(pngFilePath, "PNG")


    def decodeArray(self,npyFilePath):
        array = numpy.load(npyFilePath)

        width,height,depth = array.shape

        text = list()

        for w in range(width):
            row = list()
            for h in range(height):
                channels = array[w,h,:-1]
                index = int(numpy.sum(channels*self.decodeIndex))
                row.append(self.decodeToTextMap[index])

            text.append(''.join(row))

        return text


    def decodeRGB(self,imageFilePath):
        '''
        Read a RGB and convert to a text alignment
        :param filename:
        :return:
        '''

        image = Image.open(imageFilePath)

        pixels = numpy.array(image.getdata())
        text = list()

        width,height = image.size
        depth = 4

        pixels = numpy.reshape(pixels,(height,width,depth))

        for h in range(height):
            row = list()
            for w in range(width):
                r,g,b,a = self.RGBtoBinary(pixels[h][w])
                row.append(self.RGBtoText[r][g][b])

            text.append(''.join(row))

        return text


decoder = Decoder("/Users/kishwar/Kishwar/pileup-output/run-01122017195751/chr3/summary_3_100921-113629.csv",
                  localPath="/Users/kishwar/Kishwar/software/pileupGenerator/src/output")


decoder.convertArraysToPNG()
