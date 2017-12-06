import numpy
from PIL import Image
import os

class Decoder:
    def __init__(self,summaryFilename,localPath=None):
        self.summaryFile = summaryFilename
        self.localPath = localPath

        # matches
        self.decodeToSNPMap = ['M', 'A', 'C', 'G', 'T', 'MM', 'I', 'D', '_']
        self.decodeToTextMap = ['.', 'A', 'C', 'G', 'T', 'MM', '-', 'D', '_']
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
               self.noneChar: [0,  0,  0]}


    def reshapeArray(self,inArray,outShape):

        outArray = numpy.reshape(inArray,outShape)

        return outArray

    #
    # def decodeArraysToStdout(self):
    #     with open(self.summaryFile,'r') as file:
    #         for line in file:
    #             arrayPath, label, h, w, d = line.split(',')
    #
    #             h = int(h)
    #             w = int(w)
    #             d = int(d)
    #
    #             if self.localPath is not None:
    #                 arrayFilename = arrayPath.split('/')[-1]
    #                 arrayPath = self.localPath + arrayFilename
    #
    #             if os.path.exists(arrayPath):
    #                 pileupText = self.decodeArray(arrayPath)
    #
    #                 print(arrayPath.split('/')[-1])
    #                 print(label)
    #                 for r, row in enumerate(pileupText):
    #                     print(row)
    #             else:
    #                 print("WARNING file not found: ",arrayPath)
    #                 pass


    def convertFlattenedPNGsToPNG(self):
        with open(self.summaryFile,'r') as file:
            for line in file:
                # print(line)

                pngPath, label, h, w, d = line.split(',')

                h = int(h)
                w = int(w)
                d = int(d)

                if self.localPath is not None:
                    pngFilename = pngPath.split('/')[-1]
                    pngPath = self.localPath + pngFilename

                if os.path.exists(pngPath):
                    self.convertFlattenedPNGToRGBA(pngPath,h,w,d)
                else:
                    print("WARNING file not found: ",pngPath)
                    pass


    # def convertArraysToPNG(self):
    #     with open(self.summaryFile,'r') as file:
    #         for line in file:
    #             arrayPath, label, h, w, d = line.split(',')
    #
    #             if self.localPath is not None:
    #                 arrayFilename = arrayPath.split('/')[-1]
    #                 arrayPath = self.localPath + arrayFilename
    #
    #             if os.path.exists(arrayPath):
    #                 self.convertArrayToRGBA(arrayPath)
    #             else:
    #                 print("WARNING file not found: ",arrayPath)
    #                 pass


    # def decodeRGBstoStdout(self):
    #     with open(self.summaryFile,'r') as file:
    #         for line in file:
    #             imagePath, label = line.split(',')
    #
    #             if self.localPath is not None:
    #                 imageFilename = imagePath.split('/')[-1]
    #                 imagePath = self.localPath + imageFilename
    #
    #             if os.path.exists(imagePath):
    #                 pileupText = self.decodeRGB(imagePath)
    #
    #                 print(imagePath.split('/')[-1])
    #                 for r, row in enumerate(pileupText):
    #                     print(label[r], row)
    #             else:
    #                 print("WARNING file not found: ",imagePath)
    #                 pass


    def RGBtoBinary(self,rgb):
        return [int(value/255) for value in rgb]


    def RGBtoSortingKey(self,rgb):
        i1,i2,i3 = self.RGBtoBinary(rgb)
        code = self.RGBtoText[i1][i2][i3]
        return self.sortingKey[code]


    def convertFlattenedPNGToRGBA(self,inputPNGFilePath,height,width,depth):
        array = numpy.asarray(Image.open(inputPNGFilePath))

        # array = numpy.array(image.getdata()) #.reshape(h,w*d)

        # print(inputPNGFilePath)
        # print('h',height,'w',width,'d',depth)
        # print("raw",array.shape)
        array = self.reshapeArray(array,(height,width,depth))
        # print(array.shape)

        image = Image.new("RGBA", (height, width))
        pixels = image.load()

        for w in range(width):
            for h in range(height):
                channels = array[h, w, :-1]
                quality = array[h, w, -1]  # q is always last channel

                # print(h,w,inputPNGFilePath,channels)
                # print("q",quality)
                # print(numpy.where(channels==255)[0])

                decodeIndex = numpy.where(channels==255)[0][0]

                # print(decodeIndex)

                snp = self.decodeToSNPMap[decodeIndex]

                rgb = self.SNPtoRGB[snp]
                rgba = rgb+[quality]
                # print(quality)
                rgba = tuple(map(int, rgba))
                pixels[h, w] = rgba

        pngFilePath = inputPNGFilePath.split('.')[0]+"_decoded.png"
        image.save(pngFilePath, "PNG")


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

        pngFilePath = npyFilePath.split('.')[0] + ".png"
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



# decoder = Decoder("/Users/saureous/code/Git/pileupGenerator/data/labels/summary_1:142535300-142535600.csv",
#                   localPath="/Users/saureous/code/Git/pileupGenerator/data/labels/")

# decoder = Decoder("/Users/saureous/code/Git/pileupGenerator/data/coverage/summary_3:100000-200000.csv",
#                   localPath="/Users/saureous/code/Git/pileupGenerator/data/coverage/")

decoder = Decoder("/Users/saureous/code/Git/pileupGenerator/data/errors/summary_3:100000-200000.csv",
                  "/Users/saureous/code/Git/pileupGenerator/data/errors/")

decoder.convertFlattenedPNGsToPNG()
