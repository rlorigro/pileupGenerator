from pysam import VariantFile
from collections import defaultdict,Counter
import sys

confusion_sites_filename = "/Users/saureous/data/decoded_mismatch-short.txt"
vcf_filename = "/Users/saureous/data/NA12878_S1.genome.vcf.gz"
vcf_confident_filename = "/Users/saureous/data/NA12878_S1_confident.genome.vcf.gz"

# vcf_confident = VariantFile(vcf_confident_filename)
# vcf_all = VariantFile(vcf_filename)


class ConfusionValidator:
    def __init__(self,vcf_all_filename,vcf_confident_filename,confusion_sites_filename):
        self.vcf_all = VariantFile(vcf_all_filename)
        self.vcf_confident = VariantFile(vcf_confident_filename)
        self.confusion_sites_filename = confusion_sites_filename

        self.confusion_data = list()

        self.default_data_keys = ["region",
                                  "position",
                                  "genotype_label",
                                  "genotype_predict",
                                  "genotype_vcf",
                                  "reference_char",
                                  "reference_char_vcf",
                                  "vcf_label_mismatch_all",
                                  "vcf_reference_mismatch_all",
                                  "vcf_quality",
                                  "in_confident",
                                  "in_all",
                                  "is_insert",
                                  "frequencies",
                                  "n_characters",
                                  "n_entries_all",
                                  "record_text_all",
                                  "record_text_confident"]

        self.default_key_index = {key:self.default_data_keys.index(key) for key in self.default_data_keys}

        # self.default_data_values = [None]*len(self.default_data_keys)

        self.genotype_class_labels = {"hom_ref":1,
                                      "het_alt":2,
                                      "hom_alt":3,
                                      "none":0}

        self.insert_char = '_'


    def collect_data(self):
        '''
        Open the confusion data file and store all entries as a 3d dict.

        e.g.: data[chr1][12345][genotype_label]
                               [genotype_predict]
                               [in_confident]
                               [in_all]
                               [is_insert]
                               [frequencies]
                               [n_characters]
        :return:
        '''

        is_data = False
        region = None
        position = None
        with open(self.confusion_sites_filename,'r') as file:
            for l,line in enumerate(file):
                if is_data:
                    data = [None]*len(self.default_data_keys)

                    genotype_label,genotype_predict,a,column = line.strip().split()
                    reference_char = column[0]
                    frequencies = list(Counter(column[1:]).items())

                    in_all, \
                    n_entries_all, \
                    vcf_label_mismatch_all, \
                    genotype_vcf_all, \
                    vcf_reference_mismatch_all, \
                    reference_char_vcf_all, \
                    record_text_all, \
                    vcf_quality_all = self.check_vcf_for_entry(region=region,
                                                               position=position,
                                                               reference_char=reference_char,
                                                               genotype_label=genotype_label,
                                                               confident=False)

                    in_confident, \
                    n_entries_confident, \
                    vcf_label_mismatch_confident, \
                    genotype_vcf_confident, \
                    vcf_reference_mismatch_confident, \
                    reference_char_vcf_confident, \
                    record_text_confident,\
                    vcf_quality_confident = self.check_vcf_for_entry(region=region,
                                                                     position=position,
                                                                     reference_char=reference_char,
                                                                     genotype_label=genotype_label,
                                                                     confident=True)

                    if reference_char == '_':
                        is_insert = True
                    else:
                        is_insert = False

                    # print(genotype_label,genotype_predict,is_insert,reference_char,frequencies)

                    data[self.default_key_index["region"]] = region
                    data[self.default_key_index["position"]] = position
                    data[self.default_key_index["genotype_vcf"]] = genotype_vcf_all
                    data[self.default_key_index["genotype_label"]] = genotype_label
                    data[self.default_key_index["genotype_predict"]] = genotype_predict
                    data[self.default_key_index["reference_char"]] = reference_char
                    data[self.default_key_index["reference_char_vcf"]] = reference_char_vcf_all
                    data[self.default_key_index["is_insert"]] = is_insert
                    data[self.default_key_index["in_all"]] = in_all
                    data[self.default_key_index["in_confident"]] = in_confident
                    data[self.default_key_index["n_entries_all"]] = n_entries_all
                    data[self.default_key_index["vcf_label_mismatch_all"]] = vcf_label_mismatch_all
                    data[self.default_key_index["vcf_reference_mismatch_all"]] = vcf_reference_mismatch_all
                    data[self.default_key_index["record_text_all"]] = record_text_all
                    data[self.default_key_index["record_text_confident"]] = record_text_confident
                    data[self.default_key_index["vcf_quality"]] = vcf_quality_all
                    data[self.default_key_index["n_characters"]] = len(frequencies)

                    # print(data)

                    self.confusion_data.append(data)

                    is_data = False     # switch off


                elif line.startswith("chr"):
                    region,position = line.strip().split(':')
                    region = int(region[3])
                    position = int(position)

                    is_data = True      # switch on

                else:
                    pass    # blank line

                if l%1000 == 0:
                    sys.stderr.write("Completed %d lines\n"%l)


    def add_data(self, region, position, key, value):
        '''
        Store data relevant to a confusion site for a given data key (from class template) and value.
        :param region:
        :param position:
        :param key:
        :param value:
        :return:
        '''

        if region not in self.confusion_data:
            self.confusion_data[region] = dict()

        if position not in self.confusion_data:
            self.confusion_data[region][position] = {key:value for key,value in zip(self.default_data_keys,self.default_data_values)}

        if key not in self.confusion_data[region][position]:
            exit("ERROR: invalid key detected: %s"%(key))
        else:
            self.confusion_data[region][position][key] = value


    def check_vcf_for_entry(self,region,position,reference_char,genotype_label,confident):
        if confident:
            vcf = self.vcf_confident
        else:
            vcf = self.vcf_all

        record = None
        record_text = ''
        vcf_label_mismatch = False
        vcf_reference_mismatch = False
        vcf_found = False
        reference_char_vcf = ''
        genotype_vcf = -1
        vcf_quality = -1
        n_entries = 0

        contig = "chr"+str(region)
        entries = vcf.fetch(contig=contig,start=position-1,stop=position)
        for record in entries:
            genotype_vcf = int(self.get_class_for_genotype(self.get_genotype(record)))
            # print(genotype_vcf)
            if genotype_vcf != 1:   # ignore 0/0 genotypes, because why are they even in the VCF in the first place??
                reference_char_vcf = record.ref
                # print(region,position,record.ref)
                record_text += str(record).strip() + '\t'
                vcf_quality = record.qual if record.qual is not None else -1
                n_entries += 1

        record_text = record_text.replace(',',';').replace('\t',"    ")

        # print(record_text)

        if int(genotype_label) != genotype_vcf and genotype_vcf != 1:
            vcf_label_mismatch = True
            # print("WARNING: incorrect labelling for %s, %s"%(contig,position))
            # print('\t',record)

        if reference_char != self.insert_char \
                and (reference_char not in reference_char_vcf or reference_char != reference_char_vcf) \
                and genotype_vcf != 1:

            vcf_reference_mismatch = True
            # print("WARNING: pileup reference allele (%s) not found in VCF (%s) for %s, %s"%(reference_char,reference_char_vcf,region,position))
            # print('\t',record)

        # n_entries = len(list(entries))
        if n_entries == 0:
            vcf_found = False
            # print("WARNING: No VCF entry found for %s, %s"%(contig,position))
        else:
            vcf_found = True
            # if n_entries > 1:
                # print("WARNING: Multiple entries found for %s, %s"%(contig,position))

        return vcf_found, n_entries, vcf_label_mismatch, genotype_vcf, vcf_reference_mismatch, reference_char_vcf, record_text, vcf_quality


    def get_genotype(self,vcf_record_obj):
        return str(vcf_record_obj).rstrip().split('\t')[-1].split(':')[0].replace('/', '|').replace('\\', '|').split('|')


    def get_class_for_genotype(self, genotype_field):
        if genotype_field[0] == '.':
            return self.genotype_class_labels["hom_ref"]
        if genotype_field[0] == genotype_field[-1]:
            if genotype_field[0] == '0':
                return self.genotype_class_labels["hom_ref"]  # homozygous reference
            else:
                return self.genotype_class_labels["hom_alt"]  # homozygous alt
        else:
            return self.genotype_class_labels["het_alt"]  # heterozygous (single or double alt)


    def save_data(self,output_filename,sort_key=None,cutoff=sys.maxsize,region=1):
        if sort_key is not None:
            data = sorted(self.confusion_data, key=lambda x: x[self.default_key_index[sort_key]])
        else:
            data = self.confusion_data

        # print(data[:5])

        with open(output_filename,'w') as output_file:
            output_file.write(','.join(self.default_data_keys)+'\n')
            for e,entry in enumerate(data):
                output_file.write(','.join(list(map(str,entry)))+'\n')

                if e == cutoff:
                    break


validator = ConfusionValidator(vcf_all_filename=vcf_filename,
                               vcf_confident_filename=vcf_confident_filename,
                               confusion_sites_filename=confusion_sites_filename)

validator.collect_data()
validator.save_data("/Users/saureous/data/confusion/confusion_data.csv",cutoff=5000)