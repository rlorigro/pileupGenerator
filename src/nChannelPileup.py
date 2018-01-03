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

    def __init__(self, sam, fasta, chromosome, query_start, variant_position, flank_length, output_filename, coverage_cutoff, map_quality_cutoff, window_cutoff, coverage_threshold):
        self.length = flank_length*2+1
        self.flank_length = flank_length
        self.coverage_cutoff = coverage_cutoff
        self.output_filename = output_filename
        self.variant_position = variant_position
        self.query_start = query_start
        self.query_end = query_start+self.length-1
        self.chromosome = chromosome
        self.map_quality_cutoff = map_quality_cutoff
        self.window_cutoff = window_cutoff
        self.coverage_threshold = coverage_threshold

        # pysam fetch reads
        self.local_reads = sam.fetch("chr"+self.chromosome, start=self.query_start, end=self.query_end)

        # pyfaidx fetch reference sequence for query region
        self.ref_sequence = fasta.get_seq(name="chr"+self.chromosome, start=self.query_start, end=self.query_end)
        # self.outputRefSequence = copy.deepcopy(self.refSequence)
        self.reference_encoding = list()

        self.delta_ref  = [1,0,1,0,0,0,0,0,0,0,0]    # key for whether reference advances
        self.delta_read = [1,1,0,0,0,0,0,0,0,0,0]    # key for whether read sequence advances
                                                    #  ['M','I','D','N','S','H','P','=','X','B','NM']

        self.insert_char = '-'
        self.none_char = '_'       # character to use for empty positions in the text pileup

        self.code_template = {'M': [1.0,0.0,0.0,0.0,0.0,0.0,0.0],  # Match channel, 1=mismatch 0=match
                              'A': [0.0,1.0,0.0,0.0,0.0,0.0,0.0],
                              'C': [0.0,0.0,1.0,0.0,0.0,0.0,0.0],
                              'G': [0.0,0.0,0.0,1.0,0.0,0.0,0.0],
                              'T': [0.0,0.0,0.0,0.0,1.0,0.0,0.0],
                              'I': [0.0,0.0,0.0,0.0,0.0,1.0,0.0],
                              'D': [0.0,0.0,0.0,0.0,0.0,0.0,1.0],
                              'N': [0.0,0.0,0.0,0.0,0.0,0.0,1.0],  # redundant in case of read containing 'N'... should this be independent?
                   self.none_char: [0.0,0.0,0.0,0.0,0.0,0.0,0.0]}

        self.empty_channel_vector = [0.0 for i in range(len(self.code_template['M']))]

        self.channel_key = {'M':0,
                            'A':1,
                            'C':2,
                            'G':3,
                            'T':4,
                            'I':5,
                            'D':6,
                            'N':6,
                            self.none_char:None}

        self.decode_map = ['M', 'A', 'C', 'G', 'T', self.insert_char, 'D', self.none_char]
        self.decode_index = list(range(len(self.decode_map)))

        self.none_alpha = [1]
        self.insert_alpha = [1]
        self.reference_alpha = [1]

        self.insert_columns = dict()
        self.pileup_image = [[self.code_template[self.none_char]+self.none_alpha]*(self.window_cutoff) for j in range(self.coverage_cutoff)]  # mixing tuples (below) and lists... bad!
        self.cigar_legend = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B', '?']  # ?=NM

        self.cigar_tuples = None
        self.read_sequence = None
        self.ref_start_absolute = None
        self.ref_end_absolute = None
        self.ref_alignment_positions = None
        self.ref_index = None
        self.read_index = None
        self.read_index_relative = None
        self.ref_index_relative = None
        self.pileup_row_ends = dict()
        self.pack_row_mapping = dict()
        self.ref_anchor_positions = list()
        self.coverage_per_column = defaultdict(int)
        self.array_shape = None
        self.n_inserts_before_variant_site = 0


    def generate_decode_map(self):
        '''
        Generate inverse of code_template automatically, from any edited version of the code_template dictionary.
        Assumes only a single 1.0 exists per code... :(
        :return:
        '''

        decode_map = [None for n in range(len(self.code_template.keys()))]

        for key in self.code_template:
            vector = self.code_template[key]
            i = vector.index(1)

            decode_map[i] = key


    def generate_read_mapping(self, r, startIndex):
        '''
        Find the topmost available space to insert a read into the pileup
        :param r:
        :param startIndex:
        :return:
        '''

        i = 0
        unmapped = True
        for item in sorted(self.pileup_row_ends.items(), key=lambda x: x[1]):
            i += 1
            end = item[1]

            if end < startIndex:
                self.pack_row_mapping[r] = item[0]   # store for use during insert homogenization
                unmapped = False
                break

        if unmapped:
            self.pack_row_mapping[r] = i


    def get_match_encoding(self):
        nt = self.read_sequence[self.read_index]            # find read nucleotide
        ntRef = self.ref_sequence[self.ref_index_relative]  # find ref nucleotide

        if nt != ntRef and nt != 'N':                       # if mismatch AND NOT A SEQUENCING ERROR
            encoding = list(self.code_template['M'])        # - encode as a mismatch
        else:                                               # else
            encoding = list(self.empty_channel_vector)      # - encode as an empty template

        nt_channel = self.channel_key[nt]
        encoding[nt_channel] = 1.0                          # encode the nucleotide identity
        return encoding


    def get_delete_encoding(self):
        encoding = list(self.code_template['D'])    # encode as delete
        encoding[self.channel_key['M']] = 1.0       # set mismatch channel to 1
        return encoding


    def get_refskip_encoding(self):
        encoding = list(self.code_template['N'])  # encode as N, which doesn't have its own channel
        return encoding


    def add_pileup_entry(self, r, snp, quality):
        '''
        For a given read and SNP, add the corresponding encoding to the pileup array
        :param r:
        :param snp:
        '''

        if snp < 4 and snp != 1:
            encoding = None
            coverage = None
            self.pileup_row_ends[self.pack_row_mapping[r]] = self.ref_index_relative  # update last entry in pileupEnds

            if snp == 0:                                    # ----- MATCH --------
                encoding = self.get_match_encoding()
                coverage = 1

            elif snp == 2:                                  # ----- DELETE -------
                encoding = self.get_delete_encoding()
                coverage = 1

            elif snp == 3:                                  # ----- REFSKIP ------
                encoding = self.get_refskip_encoding()
                coverage = 0

            self.coverage_per_column[self.ref_index_relative] += coverage

            encoding = encoding+[quality]  # append the quality Alpha value and store as tuple
            self.pileup_image[self.pack_row_mapping[r]][self.ref_index_relative] = encoding      # Finally add the code to the pileup


    def add_insert_entry(self, readCharacters, n, qualities, r):
        '''
        Given a column index (i), create a dictionary of (read index):characters
        :param i:
        :param readCharacters:
        :param n:
        :param qualities:
        :return:
        '''
        i = self.ref_index_relative
        if i not in self.insert_columns:
            self.insert_columns[i] = []    # insert entries are a list of n columns where n is the length of the longest insert

        for c, character in enumerate(readCharacters):
            if c >= len(self.insert_columns[i]):   # no insert column exists yet
                self.insert_columns[i] = self.insert_columns[i]+[[self.code_template['I']+self.insert_alpha]*(self.coverage_cutoff+1)]

            encoding = list(self.code_template[character])

            encoding[self.channel_key['M']] = 1.0        # set mismatch channel to 1
            encoding[self.channel_key['I']] = 1.0        # set insert channel to 1

            self.insert_columns[i][c][self.pack_row_mapping[r]+1] = encoding+[qualities[c]]    # store quality-appended vector


    def get_pileup_encoding(self, cigarCode, refCharacter, readCharacter):
        if cigarCode == 'D':
            return 'D'
        elif readCharacter == refCharacter:
            return 'M'
        else:
            return readCharacter


    def iterate_reads(self):
        '''
        For all the reads mapped to the region specified, (subsample,) collect their alignment data, and build a draft
        pileup
        '''

        pileup_iterator_index = 0

        for r, read in enumerate(self.local_reads):
            if read.mapping_quality < self.map_quality_cutoff:
                continue
            self.parse_read(pileup_iterator_index, read)
            pileup_iterator_index += 1


    def calculate_quality(self, base_quality, map_quality):
        '''
        Approximate the overall quality of the given read character based on its mapping and sequence quality
        :param base_quality:
        :param map_quality:
        :return:
        '''
        quality = min(base_quality, map_quality)        # calculate minimum of quality for base and read quality
        quality = (1-(10 ** ((quality)/-10)))           # convert to probability and scale to 0-255 for pixel encoding

        return quality


    def parse_read(self, r, read):
        '''
        Create a draft pileup representation of a read using cigar string and read data. This method iterates through
        a read and its cigar string from the start (so it needs to be optimized for long reads).
        :param r:
        :param read:
        :return:
        '''

        ref_positions = read.get_reference_positions()
        read_qualities = read.query_qualities
        map_quality = read.mapping_quality

        first_pass = True

        if len(ref_positions) == 0:
            # sys.stderr.write("WARNING: read contains no reference alignments: %s\n" % read.query_name)
            pass
        else:
            self.ref_start_absolute = ref_positions[0] + 1
            self.cigar_tuples = read.cigartuples

            ref_positions = None  # dump

            self.read_sequence = read.query_alignment_sequence
            self.ref_index = 0              # absolute reference position that the read aligns to
            self.read_index = 0             # index of the queried read sequence
            self.ref_index_relative = 0     # index of the reference sequence to compare mismatches with

            if self.ref_start_absolute >= self.query_start:                         # if the read starts during the query region
                self.ref_index_relative = self.ref_start_absolute-self.query_start  # add the difference to ref index

            for c, entry in enumerate(self.cigar_tuples):
                snp = entry[0]
                n = entry[1]

                self.absolute_position = self.ref_index+self.ref_start_absolute
                if self.absolute_position < self.query_start-n:           # skip by blocks if not in query region yet
                    self.ref_index += self.delta_ref[snp]*n
                    self.read_index += self.delta_read[snp]*n

                elif self.absolute_position >= self.query_start-n and self.absolute_position < self.query_end:    # if within 1 cigar block of the query start, iterate 1 by 1
                    if snp == 1:    # insert
                        if self.absolute_position > self.query_start:   # if within query region
                            read_characters = self.read_sequence[self.read_index:self.read_index+n]
                            base_qualities = read_qualities[self.read_index:self.read_index+n]

                            qualities = [self.calculate_quality(base_quality, map_quality) for base_quality in base_qualities]

                            self.add_insert_entry(read_characters, n, qualities, r)

                        self.read_index += self.delta_read[snp]*n


                    elif snp < 4:   # not an insert
                        for i in range(n):

                            if self.absolute_position >= self.query_start and self.absolute_position <= self.query_end:  # if within query region
                                if first_pass:
                                    self.generate_read_mapping(r, self.ref_index_relative)
                                    first_pass = False

                                    if self.pack_row_mapping[r] >= self.coverage_cutoff:   # skip unmapped read if exceeds coverage cutoff parameter
                                        return

                                base_quality = read_qualities[self.read_index]

                                quality = self.calculate_quality(base_quality, map_quality)

                                self.add_pileup_entry(r, snp, quality)

                                self.ref_index_relative += self.delta_ref[snp]
                                self.ref_index += self.delta_ref[snp]
                                self.read_index += self.delta_read[snp]
                            else:
                                self.ref_index += self.delta_ref[snp]
                                self.read_index += self.delta_read[snp]

                            self.absolute_position = self.ref_index+self.ref_start_absolute

                else:  # stop iterating after query region
                    break


    def encode_reference(self):
        '''
        Encode the reference sequence into RGB triplets and add add it as the header line to the pileup RGB
        :return:
        '''

        for c,character in enumerate(self.ref_sequence):
            encoding = self.code_template[character]+self.reference_alpha
            self.reference_encoding.append(encoding)


    def clean_insert_columns(self, i):
        '''
        If an insert column position doesn't have read data to its left, it should be a None,
        :param i:
        :return:
        '''

        for c in range(len(self.insert_columns[i])):
            for r in range(self.coverage_cutoff+1):
                prev_entry = self.pileup_image[r][i-1][:-1]   # exclude quality

                try:
                    index = prev_entry.index(1)
                except(ValueError):
                    index = len(self.decode_map)-1

                prev_entry = self.decode_map[index]

                if prev_entry == self.none_char:    # if previous entry decodes to None character
                    self.insert_columns[i][c][r] = self.code_template[self.none_char]+self.none_alpha  # this entry should be None


    def save_pileup_array(self, filename):
        '''
        Save the pileup binary array as a bitmap image using a gray color map
        :param filename:
        :return:
        '''

        if not filename.endswith(".png"):
            filename += ".png"

        self.pileup_image = self.pileup_image[:self.coverage_cutoff]
        self.pileup_image = [row[:self.window_cutoff] for row in self.pileup_image]

        self.encode_reference()
        self.pileup_image = [self.reference_encoding]+self.pileup_image

        pileup_array = [[None for m in range(self.window_cutoff)] for n in range(self.coverage_cutoff)]

        # count insert columns that appear before variant column
        for i in range(0, self.flank_length):
            if i in self.insert_columns:
                self.n_inserts_before_variant_site += len(self.insert_columns[i])

        image_iterator = -self.n_inserts_before_variant_site
        i = 0
        while image_iterator < self.window_cutoff:
            if i in self.insert_columns:
                insert_rows = self.insert_columns[i]
                self.clean_insert_columns(i)

                for insert in insert_rows:
                    for j in range(self.coverage_cutoff):
                        if image_iterator >= 0:     # offset to ensure variant is centered
                            pileup_array[j][image_iterator] = insert[j] if j < len(insert) else self.code_template[self.none_char]+self.none_alpha

                    image_iterator += 1
                    self.ref_anchor_positions.append(i-1)
                    if image_iterator >= self.window_cutoff:
                        break

            if image_iterator >= self.window_cutoff:
                break

            for j in range(self.coverage_cutoff):
                if image_iterator >= 0:             # offset to ensure variant is centered
                    if j < self.coverage_cutoff:
                        pileup_array[j][image_iterator] = self.pileup_image[j][i] if i < len(self.pileup_image[j]) else self.code_template[self.none_char]+self.none_alpha
                    else:
                        pileup_array[j][image_iterator] = self.code_template[self.none_char]+self.none_alpha

            if i < len(self.ref_sequence):
                self.ref_anchor_positions.append(i)
            else:
                self.ref_anchor_positions.append(len(self.ref_sequence)-1)

            image_iterator += 1
            i += 1


        pileup_array = numpy.array(pileup_array)

        # print(pileup_array.shape)

        #------
        pileup_array_2d = pileup_array.reshape((pileup_array.shape[0],-1))  # unfold the stack of channels into 2d, so it can be compressed as a PNG
        self.array_shape = pileup_array.shape

        misc.imsave(self.output_filename+".png", pileup_array_2d, format="PNG")
        #------

        #------
        # pileup_array_2d = pileup_array.reshape((pileup_array.shape[0],-1))
        # sparseMatrix = sparse.csc_matrix(pileup_array_2d)
        #
        # arrayFilename = self.outputFilename + ".npz"
        # sparse.save_npz(arrayFilename,sparseMatrix)
        #------

        # numpy.save(self.outputFilename,pileup_array) # <- waste of space


class PileupGenerator:
    '''
    Creates pileups of aligned reads given a SAM/BAM alignment and FASTA reference
    '''

    def __init__(self, alignment_file, reference_file):    #, smry_ref_pos_file_writer):
        self.sam = pysam.AlignmentFile(alignment_file, "rb")
        self.fasta = Fasta(reference_file, as_raw=True, sequence_always_upper=True)

    def generate_pileup(self, chromosome,
                        position,
                        flank_length,
                        output_filename,
                        coverage_cutoff,
                        map_quality_cutoff,
                        window_cutoff,
                        coverage_threshold):
        '''
        Generate a pileup centered on a given position
        :param chromosome:
        :param position:
        :param flank_length:
        :param output_filename:
        :param coverage_cutoff:
        :param map_quality_cutoff:
        :param window_cutoff:
        :param coverage_threshold:
        :return:
        '''
        query_start = position-flank_length+1

        chromosome = str(chromosome)


        startTime = datetime.now()
        pileup = Pileup(self.sam,
                        self.fasta,
                        chromosome=chromosome,
                        query_start=query_start,
                        variant_position=position+1,
                        flank_length=flank_length,
                        output_filename=output_filename,
                        window_cutoff=window_cutoff,
                        coverage_cutoff=coverage_cutoff,
                        map_quality_cutoff=map_quality_cutoff,
                        coverage_threshold=coverage_threshold)


        # print(datetime.now() - startTime, "initialized")
        pileup.iterate_reads()
        # pileup.generateDecodeMap()
        # print(datetime.now() - startTime, "finalized")
        pileup.save_pileup_array(output_filename)
        # print(datetime.now() - startTime, "encoded and saved")
        # print()

        return pileup.query_start, pileup.ref_anchor_positions, pileup.array_shape


