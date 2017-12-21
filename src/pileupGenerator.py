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
"""
This program takes an alignment file (bam) and a reference file
to create a sparse bitmap representation of the pileup. It uses
the SamPileupBMP class and encodes each base in pileup to 6 binary
bits. It creates a large binary sparse matrix too.
"""

allVariantRecord = {}
# subregion = ':662800-663000'
# subregion = ':725000-725300'
# subregion = ':1-200000'
# subregion = ':180310000-180340000'
subregion = ''
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
    # dictionaryKeys = ["genotypeClass","InsertLength","isDelete","isMismatch","uncorrectedGenotypeClass"]

    for rec in vcf_in.fetch(region="chr"+vcf_region+subregion):
        gtField = getGTField(rec)   # genotype according to the vcf
        genotypeClass = getClassForGenotype(gtField)
        uncorrectedGenotypeClass = None

        if genotypeClass != 1 and rec.qual is not None and rec.qual > qualityCutoff:
            alleles = rec.alleles

            insertLength = None
            deleteLength = None
            isMismatch = False
            isDelete = False

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

            if rec.start in allVariantRecord:
                if allVariantRecord[rec.start][3] == True:  # overwrite only if preexisting entry is True
                    isMismatch = True
                if allVariantRecord[rec.start][2] == True:
                    isDelete = True

            if deleteLength is not None:
                offset = refLength - deleteLength
                for i in range(offset,offset+deleteLength):         # starting after the reference sequence ends, label as delete
                    if rec.start+i in allVariantRecord:             # if there is a preexisting record, only overwrite isDelete
                        record = allVariantRecord[rec.start+i]
                        record[2] = True
                        allVariantRecord[rec.start+i] = record
                    else:                                       # if no preexisting record, write new record
                        allVariantRecord[rec.start+i] = [genotypeClass, 0, True, False, genotypeClass]

            if longestLength == refLength or shortestLength == refLength:   # SNP (mismatch) is here
                isMismatch = True

            if not isMismatch and (insertLength is not None or deleteLength is not None):
                uncorrectedGenotypeClass = genotypeClass
                genotypeClass = 1
            else:
                uncorrectedGenotypeClass = genotypeClass

            # print(rec.start,rec.pos,isDelete,isMismatch)

            allVariantRecord[rec.start] = [genotypeClass, insertLength, isDelete, isMismatch, uncorrectedGenotypeClass]    # del will never be True at the anchor position


def getLabel(start, end):
    labelStr = ''              # draft of the labelling string, assuming no inserts
    insertLengths = dict()     # stores length of longest variant for labelling of inserts during pileup generation
    deletes = set()
    deleteGenotypes = dict()
    insertGenotypes = dict()
    mismatches = set()

    for j,i in enumerate(range(start, end)):
        if i in allVariantRecord.keys():
            gt = str(allVariantRecord[i][0])
            uncorrectedGenotypeClass = str(allVariantRecord[i][4])
            labelStr += gt

            if allVariantRecord[i][1] is not None:
                insertLengths[j] = int(allVariantRecord[i][1])
                insertGenotypes[j] = uncorrectedGenotypeClass
            if allVariantRecord[i][2] == True:
                deletes.add(j)
                deleteGenotypes[j] = uncorrectedGenotypeClass
            if allVariantRecord[i][3] == True:
                mismatches.add(j)

        else:
            labelStr += str(1)  # hom
    return labelStr, insertLengths, insertGenotypes, deletes, deleteGenotypes, mismatches


# def removeAnchorLabels(refStart,refPositions,label,mismatches,inserts,deletes):
#     print('M',mismatches)
#     print('I',deletes)
#     print('D',inserts)
#     print(label)
#
#     label = list(label)
#
#     for p,position in enumerate(refPositions):
#         # position-=1
#         # print(absolutePosition)
#         if position not in mismatches and (position in inserts or position in deletes):
#             print("TRUE",position)
#             label[position] = '1'
#
#     print(''.join(label))
#     print()
#     return ''.join(label)


def generatePileupBasedonVCF(vcf_region, vcf_subregion, bamFile, refFile, vcfFile, output_dir, window_size, window_cutoff,
                             coverage_cutoff, map_quality_cutoff, vcf_quality_cutoff, coverage_threshold, vcfFileConfident=None):
    files = list()
    cnt = 0
    start_timer = timer()
    populateRecordDictionary(vcf_region, vcfFile)

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

    p = SamPileupBMP.PileUpGenerator(bamFile, refFile)  #, smry_ref_pos_file_writer)
    prev_start = None
    prev_end = None

    if vcfFileConfident is not None:
        variant_sites = VariantFile(vcfFileConfident).fetch(region="chr"+vcf_region+vcf_subregion)
    else:
        variant_sites = VariantFile(vcfFile).fetch(region="chr"+vcf_region+vcf_subregion)

    for rec in variant_sites:
        if rec.qual is not None and rec.qual > vcf_quality_cutoff and getClassForGenotype(getGTField(rec)) != 1:
            start = rec.pos - window_size - 1
            end = rec.pos + window_size

            if prev_start is not None and prev_end is not None and rec.pos < prev_end - 10:
                continue
            else:
                prev_start = start
                prev_end = end

            labelString, insertLengths, insertGenotypes, deleteLengths, deleteGenotypes, mismatches = getLabel(start, end)

            filename = output_dir + rec.chrom + "_" + str(rec.pos)

            outputLabelString,refStartPosition,refAnchorPositions,arrayShape = p.generatePileup(chromosome=vcf_region,
                                                                                     position=rec.pos - 1,
                                                                                     flankLength=window_size,
                                                                                     coverageCutoff=coverage_cutoff,
                                                                                     windowCutoff=window_cutoff,
                                                                                     mapQualityCutoff=map_quality_cutoff,
                                                                                     outputFilename=filename,
                                                                                     label=labelString,
                                                                                     insertLengths=insertLengths,
                                                                                     insertGenotypes=insertGenotypes,
                                                                                     deleteLengths=deleteLengths,
                                                                                     deleteGenotypes=deleteGenotypes,
                                                                                     coverageThreshold=coverage_threshold,
                                                                                     mismatches=mismatches
                                                                                     )

            # print(filename)
            # print(labelString)
            # print(outputLabelString)

            # outputLabelString = removeAnchorLabels(refStart=rec.pos,
            #                                        label=outputLabelString,
            #                                        refPositions=refAnchorPositions,
            #                                        mismatches=mismatches,
            #                                        inserts=insertLengths,
            #                                        deletes=deleteLengths)


            cnt += 1
            if cnt % 1000 == 0:
                end_timer = timer()
                print(str(cnt) + " Records done", file=sys.stderr)
                print("TIME elapsed " + str(end_timer - start_timer), file=sys.stderr)

            smry.write(os.path.abspath(filename) + ".png," + str(outputLabelString) + ',' + ','.join(map(str,arrayShape))+'\n')
            row = [os.path.abspath(filename) + ".png,", refStartPosition]+refAnchorPositions
            smry_ref_pos_file_writer.writerow(row)

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
                             coverage_cutoff, map_quality_cutoff, vcf_quality_cutoff, threads, coverage_threshold,vcfFileConfident=None):
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
                                                           map_quality_cutoff, vcf_quality_cutoff,coverage_threshold,vcfFileConfident))
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
        "--vcf_confident",
        type=str,
        required=False,
        help="Confident VCF file containing SNPs and SVs. If provided, pileups will only be generated at sites from this \
             file. Non-confident sites may be found and labelled within the window of confident sites, however."
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
        default=0,
        help="Phred scaled threshold for mapping quality."
    )
    parser.add_argument(
        "--vcf_quality_cutoff",
        type=int,
        default=10,
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
                            vcfFileConfident = FLAGS.vcf_confident,
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
    #                          vcfFileConfident=FLAGS.vcf_confident,
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
