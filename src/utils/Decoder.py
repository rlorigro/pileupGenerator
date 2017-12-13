import numpy
from PIL import Image
import os

class Decoder:
    def __init__(self,summaryFilename,localPath=None):
        self.summaryFile = summaryFilename
        self.localPath = localPath

        self.RGBtoText = [[[' ', 'T'], ['G', 'D']], [['A', '_'], ['C', '.']]]

        self.sortingKey = {'M':5,
                           'A':0,
                           'C':1,
                           'G':2,
                           'T':3,
                           'I':6,
                           'D':4,
                           'N':7}


    def decodeToStdout(self):
        with open(self.summaryFile,'r') as file:
            for line in file:
                imagePath, label = line.split(',')

                if self.localPath is not None:
                    imageFilename = imagePath.split('/')[-1]
                    imagePath = self.localPath + imageFilename

                if os.path.exists(imagePath):
                    image = Image.open(imagePath)

                    pileupText = self.decodeRGB(image)

                    print(imagePath.split('/')[-1])
                    for r, row in enumerate(pileupText):
                        print(label[r], row)
                else:
                    # print("WARNING file not found: ",imagePath)
                    pass


    def RGBtoBinary(self,rgb):
        return [int(value/255) for value in rgb]


    def RGBtoSortingKey(self,rgb):
        i1,i2,i3 = self.RGBtoBinary(rgb)
        code = self.RGBtoText[i1][i2][i3]
        return self.sortingKey[code]


    def decodeRGB(self,image):
        '''
        Read a RGB and convert to a text alignment
        :param filename:
        :return:
        '''

        pixels = numpy.array(image.getdata())
        text = list()

        width,height = image.size
        depth = 4

        pixels = numpy.reshape(pixels,(height,width,depth))

        for h in range(height):
            row = ''
            for w in range(width):
                r,g,b,a = self.RGBtoBinary(pixels[h][w])

                row += self.RGBtoText[r][g][b]
            text.append(row)

        return text


# decoder = Decoder("/Users/saureous/Documents/Git/Basecaller/deePore/data/jarvis-11-19/summary-3-95086495-115041333.csv",
#                   localPath="/Users/saureous/Documents/Git/Basecaller/deePore/data/jarvis-11-19/")
#
# decoder.decodeToStdout()
