import argparse
from pysam import VariantFile
import multiprocessing
from multiprocessing import Process, Pool, TimeoutError, Array, Manager
import time
from timeit import default_timer as timer
import sys
import os
import csv
import nChannelPileup as SamPileupBMP
from vcf_handler import *
"""
This program takes an alignment file (bam) and a reference file
to create a sparse bitmap representation of the pileup. It uses
the SamPileupBMP class and encodes each base in pileup to 6 binary
bits. It creates a large binary sparse matrix too.
"""

allVariantRecord = {}
# subregion = ':1547900-1547990'
# subregion = ':9684214-9684290'
# subregion = ':63996-63997'
# subregion = ':1-200000'
subregion = ''

def getClassForGenotype(gtField):
    if gtField[0] == '.':
        return 1
    if gtField[0] == gtField[-1]:
        if gtField[0] == '0':
            return 1    # homozygous reference
        else:
            return 3    # homozygous alt
    else:
        return 2    # heterozygous (single or double alt)


# def getGTField(rec):
#     return str(rec).rstrip().split('\t')[-1].split(':')[0].replace('/', '|').replace('\\', '|').split('|')


class RecordGetter:
    def __init__(self,contig,subregion,vcf_handler):
        '''
        Take a vcf handler object and retrieve one entry at a time from the start to the end, ensuring all entries
        are visited if multiple exist at a single position. If none remain, None is returned
        :param start:
        :param stop:
        :param vcf_handler:
        '''
        self.contig = str(contig)
        self.subregion = subregion
        self.parser = vcf_handler
        self.start = None
        self.stop = None
        self.current_position = None
        self.current_records = None
        self.records = None
        self.record_positions = None

        if not self.contig.startswith("chr"):
            self.contig = "chr" + self.contig

        self.initialize_records()

    def initialize_records(self):
        self.parser.populate_dictionary(contig=self.contig, site=self.subregion)
        self.records = self.parser.get_variant_dictionary()

    def get_next(self):
        if self.current_position is None:                                       # if positions are uninitialized
            all_positions = sorted(list(self.records.keys()))
            self.record_positions = iter(all_positions)                         # initialize iterator for record positions
            self.current_position = next(self.record_positions)
            self.current_records = self.records[self.current_position]
            self.start = all_positions[0]
            self.stop = all_positions[-1]

        if self.current_position == self.stop and len(self.current_records) == 0:   # if all records at last position have been visited
            return None
        else:
            if len(self.current_records) == 0:                                  # if all records at any given position have been visited
                self.current_position = next(self.record_positions)             # find next position and corresponding records
                self.current_records = self.records[self.current_position]
            return self.current_records.pop()                                   # return last record at the position and remove from list


def generatePileupBasedonVCF(vcf_region, vcf_subregion, bamFile, refFile, vcfFile, output_dir, window_size, window_cutoff,
                             coverage_cutoff, map_quality_cutoff, vcf_quality_cutoff, coverage_threshold):
    files = list()
    cnt = 0
    start_timer = timer()
    # populateRecordDictionary(vcf_region, vcfFile, vcf_quality_cutoff)

    if vcf_subregion[0] == ':':
        vcf_subregion_filename = '_' + vcf_subregion[1:]
    else:
        vcf_subregion_filename = vcf_subregion

    smry = open(output_dir + "summary" + '_' + vcf_region + vcf_subregion_filename + ".csv", 'w')
    smry_ref_pos_file = open(output_dir + "ref_positions_" + vcf_region + vcf_subregion_filename + ".csv", 'w')
    smry_ref_pos_file_writer = csv.writer(smry_ref_pos_file)

    try:
        os.stat("tmp/")
    except:
        os.mkdir("tmp/")

    log = open("tmp/" + "log_" + vcf_region + vcf_subregion + ".txt", 'w')
    files.append(smry)
    files.append(log)
    files.append(smry_ref_pos_file)

    log.write("Reference file: \t%s\n" % refFile)
    log.write("BAM file: \t%s\n" % bamFile)
    log.write("VCF file: \t%s\n" % vcfFile)
    log.write("Region: \t%s\n" % vcf_region)
    log.write("VCF Subregion: \t%s\n" % vcf_subregion)
    log.write("Window size: \t%s\n" % window_size)
    log.write("Coverage cutoff: \t%s\n" % coverage_cutoff)
    log.write("VCF quality cutoff: \t%s\n" % vcf_quality_cutoff)
    log.write("Map quality cutoff: \t%s\n" % map_quality_cutoff)
    log.write("Coverage threshold: \t%s\n" % coverage_threshold)

    p = SamPileupBMP.PileupGenerator(bamFile, refFile)  #, smry_ref_pos_file_writer)
    prev_start = None
    prev_end = None

    vcf_handler = VCFFileProcessor(vcfFile)
    record_getter = RecordGetter(contig=vcf_region,subregion=vcf_subregion,vcf_handler=vcf_handler)

    rec = record_getter.get_next()
    while rec is not None:
        if rec.qual is not None and rec.qual > vcf_quality_cutoff and rec.type != "Hom":
            # print(rec.type)
            filename = output_dir + vcf_region + "_" + str(rec.pos)

            refStartPosition,refAnchorPositions,arrayShape = p.generate_pileup(chromosome=vcf_region,
                                                                               position=rec.pos - 1,
                                                                               flank_length=window_size,
                                                                               coverage_cutoff=coverage_cutoff,
                                                                               window_cutoff=window_cutoff,
                                                                               map_quality_cutoff=map_quality_cutoff,
                                                                               output_filename=filename,
                                                                               coverage_threshold=coverage_threshold
                                                                               )

            cnt += 1
            if cnt % 1000 == 0:
                end_timer = timer()
                print(str(cnt) + " Records done", file=sys.stderr)
                print("TIME elapsed " + str(end_timer - start_timer), file=sys.stderr)

            smry.write(os.path.abspath(filename) + ".png," + str(rec.type) + ',' + ','.join(map(str,arrayShape))+'\n')
            row = [os.path.abspath(filename) + ".png,", refStartPosition]+refAnchorPositions
            smry_ref_pos_file_writer.writerow(row)

        rec = record_getter.get_next()  # return the next record until None is returned
        # print(rec)

    for file in files:
        file.close()


def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out


def parallel_pileup_generator(vcf_region, bamFile, refFile, vcfFile, output_dir, window_size, window_cutoff,
                             coverage_cutoff, map_quality_cutoff, vcf_quality_cutoff, threads, coverage_threshold):
    all_positions = []

    vcf_handler = VCFFileProcessor(vcfFile)
    record_getter = RecordGetter(contig=vcf_region,subregion='',vcf_handler=vcf_handler)

    rec = record_getter.get_next()
    while rec is not None:
        if rec.qual is not None and rec.qual > vcf_quality_cutoff and rec.type != "Hom":
            all_positions.append(rec.pos)

    all_positions = chunkIt(all_positions, threads)
    starts = [all_positions[i][0] for i in range(len(all_positions))]
    ends = [all_positions[i][-1] for i in range(len(all_positions))]
    args = ()

    for i in range(len(starts)):
        args = args + ((starts[i], ends[i]),)
        vcf_subregion = "-"+str(starts[i])+"-"+str(ends[i])
        p = Process(target=generatePileupBasedonVCF, args=(vcf_region, vcf_subregion, bamFile, refFile, vcfFile,
                                                           output_dir, window_size, window_cutoff, coverage_cutoff,
                                                           map_quality_cutoff, vcf_quality_cutoff,coverage_threshold))
        p.start()


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="BAM file with alignments."
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="Reference corresponding to the BAM file."
    )
    parser.add_argument(
        "--vcf",
        type=str,
        required=True,
        help="VCF file containing SNPs and SVs."
    )
    parser.add_argument(
        "--region",
        type = str,
        default = "chr3",
        help="Site region. Ex: chr3"
    )
    parser.add_argument(
        "--vcf_region",
        type=str,
        default="3",
        help="Site region. Ex: chr3"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="output/",
        help="Name of output directory"
    )
    parser.add_argument(
        "--window_size",
        type=int,
        default=10,
        help="Window size of query region."
    )
    parser.add_argument(
        "--window_cutoff",
        type=int,
        default=30,
        help="Size of output image."
    )
    parser.add_argument(
        "--coverage",
        type=int,
        default=180,
        help="Read coverage, default is 50x."
    )
    parser.add_argument(
        "--coverage_cutoff",
        type=int,
        default=200,
        help="Size of output image."
    )
    parser.add_argument(
        "--map_quality_cutoff",
        type=int,
        default=10,
        help="Phred scaled threshold for mapping quality."
    )
    parser.add_argument(
        "--vcf_quality_cutoff",
        type=int,
        default=0,
        help="Phred scaled threshold for variant call quality."
    )
    parser.add_argument(
        "--coverage_threshold",
        type=int,
        default=5,
        help="Threshold below which to remove training label from pileup"
    )
    # parser.add_argument(
    #     "--parallel",
    #     type=bool,
    #     default=False,
    #     help="Option to run threaded pileup generator"
    # )
    parser.add_argument(
        "--max_threads",
        type=int,
        default=5,
        help="Number of maximum threads for this region."
    )
    FLAGS, unparsed = parser.parse_known_args()
    if FLAGS.window_size * 2 + 1 > FLAGS.window_cutoff:
        sys.stderr.write("ERROR: WINDOW CUTOFF TOO SMALL. MINIMUM SUPPORTED WINDOW SIZE: {2*WINDOW_SIZE + 1}\n")
        exit()

    # if FLAGS.parallel == True:

    parallel_pileup_generator(vcf_region=FLAGS.vcf_region,
                            bamFile = FLAGS.bam,
                            refFile = FLAGS.ref,
                            vcfFile = FLAGS.vcf,
                            output_dir = FLAGS.output_dir,
                            window_size = FLAGS.window_size,
                            window_cutoff = FLAGS.window_cutoff,
                            coverage_cutoff = FLAGS.coverage_cutoff,
                            map_quality_cutoff = FLAGS.map_quality_cutoff,
                            vcf_quality_cutoff = FLAGS.vcf_quality_cutoff,
                            coverage_threshold = FLAGS.coverage_threshold,
                            threads = FLAGS.max_threads)

    # else:
    # generatePileupBasedonVCF(vcf_region=FLAGS.vcf_region,
    #                          vcf_subregion=subregion,
    #                          bamFile=FLAGS.bam,
    #                          refFile=FLAGS.ref,
    #                          vcfFile=FLAGS.vcf,
    #                          output_dir=FLAGS.output_dir,
    #                          window_size=FLAGS.window_size,
    #                          window_cutoff=FLAGS.window_cutoff,
    #                          coverage_cutoff=FLAGS.coverage_cutoff,
    #                          map_quality_cutoff=FLAGS.map_quality_cutoff,
    #                          vcf_quality_cutoff=FLAGS.vcf_quality_cutoff,
    #                          coverage_threshold=FLAGS.coverage_threshold)


# example usage:
# python3 "src/pileupGenerator.py" --bam "data/chr3_200k.bam" --ref "data/chr3.fa" --vcf "data/NA12878_S1.genome.vcf.gz" --region "chr3" --vcf_region '3' --output_dir "data/errors/" --window_size 25 >data/out.txt
# python3 "src/pileupGenerator.py" --bam "/Users/saureous/data/chr3_200k.bam" --ref "/Users/saureous/data/chr3.fa" --vcf "/Users/saureous/data/NA12878_S1.genome.vcf.gz" --region "chr3" --vcf_region '3' --output_dir "/Users/saureous/data/test/" --window_size 25 --window_cutoff 75 >/Users/saureous/data/test/out.txt
