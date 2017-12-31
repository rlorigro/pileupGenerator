import pysam
import sys


class VCFRecord:
    def __init__(self, rec):
        self.rec_pos = rec.pos
        self.rec_qual = rec.qual
        self.rec_genotype = self._get_record_genotype(rec)
        self.rec_filter = list(rec.filter)[0]

        self.rec_not_hom = True

        if len(self.rec_genotype) == 2:
            self.rec_not_hom = not ((self.rec_genotype[0] == self.rec_genotype[1]
                                    and self.rec_genotype[0] == 0) or self.rec_filter != 'PASS')
        elif len(self.rec_genotype) == 1:
            self.rec_not_hom = not (self.rec_genotype[0] == 0)
        self.rec_chrom = rec.chrom
        self.rec_alleles = rec.alleles
        self.rec_alts = rec.alts
        self.rec_ref = rec.ref

    @staticmethod
    def _get_record_genotype(record):
        gts = [s['GT'] for s in record.samples.values()]
        return gts[0]


class VariantRecord:
    def __init__(self, pos, ref, qual, genotype_class, genotype_type, alt):
        self.pos = pos
        self.qual = qual
        self.ref = ref
        self.alt = alt
        self.type = genotype_type
        self.genotype_class = genotype_class

    def __str__(self):
        return_str = str(self.pos) + '\t' + str(int(self.qual)) + '\t' + str(self.ref) + '\t' + str(self.alt) +\
                     '\t' + str(self.type) + '\t' + str(self.genotype_class)
        return return_str


class VCFFileProcessor:
    def __init__(self, file_path):
        self.file_path = file_path
        self.vcf_records = None
        self.genotype_dictionary = {}

    def __str__(self):
        ret_str = "POS\tQUAL\tREF\tALT\tTYPE\tCLASS\t\n"
        for pos in self.genotype_dictionary.keys():
            for rec in self.genotype_dictionary[pos]:
                ret_str += str(rec) + "\n"

        return ret_str

    @staticmethod
    def get_genotype_class(rec, ref, alt):
        if len(ref) == 1 and len(alt) == 1:
            return 'SNP'
        elif len(ref) < len(alt):
            return 'IN'
        elif len(ref) > len(alt):
            return 'DEL'
        else:
            raise ValueError('INVALID GENOTYPE CLASS \n' + rec)

    @staticmethod
    def get_genotype_type(genotype):
        g_type = ""
        if len(genotype) == 2:
            if genotype[0] != genotype[1]:
                g_type = "Het"
            else:
                g_type = "Hom_alt"
        elif len(genotype) == 1:
            if genotype[0] == 1:
                g_type = "Hom_alt"
        else:
            raise ValueError("INVALID GENOTYPE ENCOUNTERED" + genotype)
        return g_type

    def _initialize_dictionary(self, ref_pos):
        if ref_pos not in self.genotype_dictionary.keys():
            self.genotype_dictionary[ref_pos] = []

    def _update_dictionary(self, variant_record):
        self._initialize_dictionary(variant_record.pos)
        self.genotype_dictionary[variant_record.pos].append(variant_record)

    def _process_genotype_by_class(self, rec, genotype_class, genotype_type, alt):
        if genotype_class == 'SNP':
            variant_record = VariantRecord(rec.rec_pos, rec.rec_ref, rec.rec_qual, genotype_class, genotype_type, alt)
            self._update_dictionary(variant_record)
        elif genotype_class == 'DEL':
            for pos in range(0, len(rec.rec_ref), 1):
                if pos >= len(alt):
                    variant_record = VariantRecord(rec.rec_pos + pos, rec.rec_ref[pos], rec.rec_qual, genotype_class,
                                                   genotype_type, alt='*')
                    self._update_dictionary(variant_record)
        elif genotype_class == 'IN':
            # for pos in range(0, len(rec.rec_ref), 1):
                # if rec.rec_ref[pos] != alt[pos]:
                    # raise ValueError("WEIRD INSERT", rec.rec_pos, rec.rec_ref, alt)
            pos = rec.rec_pos + len(rec.rec_ref) - 1
            ref = rec.rec_ref[-1]
            altr = alt[len(rec.rec_ref):]
            variant_record = VariantRecord(pos, ref, rec.rec_qual, genotype_class, genotype_type, altr)
            self._update_dictionary(variant_record)

    def _parse_through_records(self, vcf_record):
        for alt in vcf_record.rec_alts:
            genotype_class = self.get_genotype_class(vcf_record, vcf_record.rec_ref, alt)
            genotype_type = self.get_genotype_type(vcf_record.rec_genotype)
            self._process_genotype_by_class(vcf_record, genotype_class, genotype_type, alt)

    def _get_filtered_records(self):
        filtered_records = []
        for record in self.vcf_records:
            vcf_record = VCFRecord(record)
            if vcf_record.rec_not_hom is True and vcf_record.rec_filter == 'PASS':
                filtered_records.append(vcf_record)
        return filtered_records

    def _generate_dictionary_from_records(self, records):
        for record in records:
            self._parse_through_records(record)

    def get_variant_dictionary(self):
        return self.genotype_dictionary

    def test_filtered_records(self):
        vcf_out = pysam.VariantFile('-', 'w')
        filtered_records = self._get_filtered_records()

        for record in filtered_records:
            vcf_out.write(record)

        return vcf_out

    def populate_dictionary(self, contig, site):
        try:
            self.vcf_records = pysam.VariantFile(self.file_path).fetch(region=contig + site)
        except IOError:
            sys.stderr.write("VCF FILE READ ERROR")

        filtered_records = self._get_filtered_records()
        self._generate_dictionary_from_records(filtered_records)

