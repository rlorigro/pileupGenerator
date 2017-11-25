import argparse
from pysam import VariantFile
import multiprocessing
from multiprocessing import Process, Pool, TimeoutError, Array, Manager
import time
from timeit import default_timer as timer
import sys
import os
import csv
import PackedPileup as SamPileupBMP
"""
This program takes an alignment file (bam) and a reference file
to create a sparse bitmap representation of the pileup. It uses
the SamPileupBMP class and encodes each base in pileup to 6 binary
bits. It creates a large binary sparse matrix too.
"""

allVariantRecord = {}
subregion = ''
# subregion = ''
cutoffOutput = False
cutoff = 350


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


def getGTField(rec):
    return str(rec).rstrip().split('\t')[-1].split(':')[0].replace('/', '|').replace('\\', '|').split('|')


def populateRecordDictionary(vcf_region, vcfFile, qualityCutoff=60):
    vcf_in = VariantFile(vcfFile)
    for rec in vcf_in.fetch(region="chr"+vcf_region+subregion):
        gtField = getGTField(rec)   # genotype according to the vcf
        genotypeClass = getClassForGenotype(gtField)

        if genotypeClass != 1 and rec.qual is not None and rec.qual > qualityCutoff:
            alleles = rec.alleles

            insertLength = None
            deleteLength = None

            altAlleles = [alleles[int(gt)] for gt in gtField if gt!='0']

            altAllelesByLength = sorted(altAlleles, key=lambda x: len(x), reverse=True)

            longestLength = len(altAllelesByLength[0])
            shortestLength = len(altAllelesByLength[-1])
            refLength = len(rec.ref)

            if '0' not in gtField:                                      # double alternate
                if genotypeClass == 2:                                  # heterozygous, must consider both alt lengths
                    if longestLength == shortestLength:
                        if longestLength > refLength:                   # insert
                            insertLength = longestLength - refLength
                        elif longestLength < refLength:
                            deleteLength = refLength - longestLength    # delete
                    else:
                        if longestLength > refLength:                   # insert
                            insertLength = longestLength - refLength
                        elif longestLength < refLength:
                            deleteLength = refLength - longestLength    # delete

                        if shortestLength > refLength:                  # also insert
                            insertLength = shortestLength - refLength
                        elif shortestLength < refLength:                # also delete
                            deleteLength = refLength - shortestLength

                else:                                                   # homozygous, only one alt
                    if longestLength > refLength:                       # insert
                        insertLength = longestLength - refLength
                    else:
                        deleteLength = refLength - longestLength

            else:                       # single alt, check length of longest allele to determine if ins or del
                if longestLength > refLength:
                    insertLength = longestLength - refLength
                elif longestLength < refLength:
                    deleteLength = refLength - longestLength

            allVariantRecord[rec.start] = (genotypeClass, insertLength, deleteLength)


def getLabel(start, end):
    labelStr = ''              # draft of the labelling string, assuming no inserts
    insertLengths = dict()     # stores length of longest variant for labelling of inserts during pileup generation
    deleteLengths = dict()

    insertCount = 0
    deleteCount = 0

    for j,i in enumerate(range(start, end)):
        if i in allVariantRecord.keys():
            labelStr += str(allVariantRecord[i][0])

            # if insertCount == 0:
            if allVariantRecord[i][1] is not None:
                insertLengths[j] = int(allVariantRecord[i][1])
            #         insertCount = allVariantRecord[i][1] - 1
            # else:
            #     insertCount -= 1
            #
            # if deleteCount == 0:
            if allVariantRecord[i][2] is not None:
                deleteLengths[j] = int(allVariantRecord[i][2])
            #         deleteCount = allVariantRecord[i][2] - 1
            # else:
            #     deleteCount -= 1

        else:
            labelStr += str(1)
    return labelStr, insertLengths, deleteLengths


def generatePileupBasedonVCF(vcf_region, vcf_subregion, bamFile, refFile, vcfFile, output_dir, window_size, window_cutoff,
                             coverage_cutoff, map_quality_cutoff, vcf_quality_cutoff):
    files = list()
    cnt = 0
    start_timer = timer()
    populateRecordDictionary(vcf_region, vcfFile)
    smry = open(output_dir + "summary" + '_' + vcf_region + vcf_subregion + ".csv", 'w')
    smry_ref_pos_file = open(output_dir + "ref_positions_" + vcf_region + vcf_subregion + ".csv", 'w')
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


    p = SamPileupBMP.PileUpGenerator(bamFile, refFile, smry_ref_pos_file_writer)

    for rec in VariantFile(vcfFile).fetch(region="chr"+vcf_region+vcf_subregion):
        if rec.qual is not None and rec.qual > vcf_quality_cutoff and getClassForGenotype(getGTField(rec)) != 1:
            start = rec.pos - window_size - 1
            end = rec.pos + window_size
            labelString, insertLengths, deleteLengths = getLabel(start, end)
            filename = output_dir + rec.chrom + "_" + str(rec.pos)
            outputLabelString = p.generatePileup(chromosome=vcf_region,
                                                 position=rec.pos - 1,
                                                 flankLength=window_size,
                                                 coverageCutoff=coverage_cutoff,
                                                 windowCutoff=window_cutoff,
                                                 mapQualityCutoff=map_quality_cutoff,
                                                 outputFilename=filename,
                                                 label=labelString,
                                                 insertLengths=insertLengths,
                                                 deleteLengths=deleteLengths
                                                 )
            cnt += 1
            if cnt % 1000 == 0:
                end_timer = timer()
                print(str(cnt) + " Records done", file=sys.stderr)
                print("TIME elapsed " + str(end_timer - start_timer), file=sys.stderr)
            smry.write(os.path.abspath(filename) + ".png," + str(outputLabelString)+'\n')

            if cutoffOutput:
                if cnt > cutoff:
                    break

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
                             coverage_cutoff, map_quality_cutoff, vcf_quality_cutoff, threads):
    all_positions = []
    for rec in VariantFile(vcfFile).fetch(region="chr"+vcf_region):
        if rec.qual is not None and rec.qual > vcf_quality_cutoff and getClassForGenotype(getGTField(rec)) != 1:
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
                                                           map_quality_cutoff, vcf_quality_cutoff,))
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
        required=False,
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
        default=100,
        help="Window size of query region."
    )
    parser.add_argument(
        "--window_cutoff",
        type=int,
        default=220,
        help="Size of output image."
    )
    parser.add_argument(
        "--coverage",
        type=int,
        default=50,
        help="Read coverage, default is 50x."
    )
    parser.add_argument(
        "--coverage_cutoff",
        type=int,
        default=50,
        help="Size of output image."
    )
    parser.add_argument(
        "--map_quality_cutoff",
        type=int,
        default=5,
        help="Phred scaled threshold for mapping quality."
    )
    parser.add_argument(
        "--vcf_quality_cutoff",
        type=int,
        default=60,
        help="Phred scaled threshold for variant call quality."
    )
    parser.add_argument(
        "--max_threads",
        type=int,
        default=10,
        help="Number of maximum threads for this region."
    )

    FLAGS, unparsed = parser.parse_known_args()
    parallel_pileup_generator(FLAGS.vcf_region,
                             FLAGS.bam,
                             FLAGS.ref,
                             FLAGS.vcf,
                             FLAGS.output_dir,
                             FLAGS.window_size,
                             FLAGS.window_cutoff,
                             FLAGS.coverage_cutoff,
                             FLAGS.map_quality_cutoff,
                             FLAGS.vcf_quality_cutoff,
                             FLAGS.max_threads)
    # generatePileupBasedonVCF(FLAGS.vcf_region,
    #                          subregion,
    #                          FLAGS.bam,
    #                          FLAGS.ref,
    #                          FLAGS.vcf,
    #                          FLAGS.output_dir,
    #                          FLAGS.window_size,
    #                          FLAGS.window_cutoff,
    #                          FLAGS.coverage_cutoff,
    #                          FLAGS.map_quality_cutoff,
    #                          FLAGS.vcf_quality_cutoff)


# example usage:
# python3 "deePore/src/utils/pileupGenerator.py" --bam "deePore/data/chr3_200k.bam" --ref "deePore/data/chr3.fa" --vcf "deePore/data/NA12878_S1.genome.vcf.gz" --region "chr3" --output_dir "deePore/data/test2/" --window_size 25 >deePore/data/stdout_final.txt
