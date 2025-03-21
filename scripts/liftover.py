import collections
import logging
import sys
from copy import copy
from io import StringIO

import vcf
from vcf.parser import _Info as VcfInfo

contig_spec = collections.namedtuple('Contig', 'id,length')
calldata_spec = collections.namedtuple('CallData', 'GT')
CHROMS = [str(i) for i in range(1, 23)] + ['X', 'Y']

# the following is an minimum vcf needed to generate new ones
# If there are no records in a vcf file, we use this template to add in mismatch sites
EXAMPLE_VCF = [
    "##fileformat=VCFv4.2",
    "##reference=file:///resources/b37/GRCh38.p12.fa",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##contig=<ID=chr2>",
    "##contig=<ID=chr6>",
    "##contig=<ID=chr15>",
    "##contig=<ID=chr19>",
    "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1",
    "chr2  21012603     .      C     T     100    PASS  .    GT  1/1",
]

def convert_hg19_vcf_to_grch37_vcf(input_vcf_file, output_vcf_file):
    """
    Liftover lifts from GRCh38 -> hg19. hg19 contigs have 'chr' prefix,
    that needs to be removed for GRCH37. Also filtering to only acceptable chrs
    autosomes, X and Y. No alternate contigs are considered
    :param input_vcf_file: hg19 vcf file liftedover from GRCh38
    :param output_vcf_file: GRCh37 contigs vcf file
    """
    chroms_to_keep = ['chr' + chrom for chrom in CHROMS]
    grch37_contigs = collections.OrderedDict()
    with open(input_vcf_file) as in_fp, open(output_vcf_file, 'w') as out_fp:
        reader = vcf.Reader(in_fp)
        # update contigs to remove 'chr' prefix
        for contig in reader.contigs:
            if contig in chroms_to_keep:
                value = contig_spec(reader.contigs[contig].id.replace('chr', ''), reader.contigs[contig].length)
                grch37_contigs.update({contig.replace('chr', ''): value})

        reader.contigs = grch37_contigs
        writer = vcf.Writer(out_fp, reader, lineterminator='\n')
        for record in reader:
            if record.CHROM in chroms_to_keep:
                # update chrom in the records
                record.CHROM = record.CHROM.replace('chr', '')
                writer.write_record(record)



AOU_MISMATCH_SITES = {
    'chr2:21012603': {
        '38_coordinates': {
            'chrom': 'chr2',
            'start': 21012603,
            'end': 21012604,
            'base': 'C',
        },
        '37_coordinates': {
            'chrom': 'chr2',
            'start': 21235475,
            'end': 21235476,
            'base': 'T',
        },
    },
    'chr6:7563750': {
        '38_coordinates': {
            'chrom': 'chr6',
            'start': 7563750,
            'end': 7563751,
            'base': 'G',
        },
        '37_coordinates': {
            'chrom': 'chr6',
            'start': 7563983,
            'end': 7563984,
            'base': 'T',
        },
    },
    'chr15:48515440': {
        '38_coordinates': {
            'chrom': 'chr15',
            'start': 48515440,
            'end': 48515441,
            'base': 'T',
        },
        '37_coordinates': {
            'chrom': 'chr15',
            'start': 48807637,
            'end': 48807638,
            'base': 'C',
        },
    },
    'chr19:55154216': {
        '38_coordinates': {
            'chrom': 'chr19',
            'start': 55154216,
            'end': 55154217,
            'base': 'C',
        },
        '37_coordinates': {
            'chrom': 'chr19',
            'start': 55665584,
            'end': 55665585,
            'base': 'A',
        },
    },
}

EPILEPSY_MISMATCH_SITES =  {
    
    'chr2:21012603':{
        '38_coordinates': {
            'chrom': 'chr2',
            'start': 21012603,
            'end': 21012604,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr2',
            'start': 21235475,
            'end': 21235476,
            'base': 'T'
        }
    },
    'chr2:27217917':{
        '38_coordinates': {
            'chrom': 'chr2',
            'start': 27217917,
            'end': 27217918,
            'base': 'T'
        },
        '37_coordinates': {
            'chrom': 'chr2',
            'start': 27440785,
            'end': 27440786,
            'base': 'C'
        }
    },
    'chr2:178756750':{
        '38_coordinates': {
            'chrom': 'chr2',
            'start': 178756750,
            'end': 178756751,
            'base': 'T'
        },
        '37_coordinates': {
            'chrom': 'chr2',
            'start': 179621477,
            'end': 179621478,
            'base': 'C'
        }
    },
    'chr6:7563750':{
        '38_coordinates': {
            'chrom': 'chr6',
            'start': 7563750,
            'end': 7563751,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr6',
            'start': 7563983,
            'end': 7563984,
            'base': 'T'
        }
    },
    'chr12:2682039':{
        '38_coordinates': {
            'chrom': 'chr12',
            'start': 2682039,
            'end': 2682040,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr12',
            'start': 2791205,
            'end': 2791206,
            'base': 'A'
        }
    },
    'chr15:48515440':{
        '38_coordinates': {
            'chrom': 'chr15',
            'start': 48515440,
            'end': 48515441,
            'base': 'T'
        },
        '37_coordinates': {
            'chrom': 'chr15',
            'start': 48807637,
            'end': 48807638,
            'base': 'C'
        }
    },
    'chr16:88773':{
        '38_coordinates': {
            'chrom': 'chr16',
            'start': 88773,
            'end': 88775,
            'base': 'AG'
        },
        '37_coordinates': {
            'chrom': 'chr16',
            'start': 138772,
            'end': 138773,
            'base': 'A'
        }
    },
    'chr17:8238389':{
        '38_coordinates': {
            'chrom': 'chr17',
            'start': 8238389,
            'end': 8238390,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr17',
            'start': 8141707,
            'end': 8141708,
            'base': 'T'
        }
    },
    'chr17:46067627':{
        '38_coordinates': {
            'chrom': 'chr17',
            'start': 46067627,
            'end': 46067628,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr17',
            'start': 44144993,
            'end': 44144994,
            'base': 'C'
        }
    },
    'chr18:55635870':{
        '38_coordinates': {
            'chrom': 'chr18',
            'start': 55635870,
            'end': 55635871,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr18',
            'start': 53303101,
            'end': 53303102,
            'base': 'C'
        }
    },
    'chr18:62138299':{
        '38_coordinates': {
            'chrom': 'chr18',
            'start': 62138299,
            'end': 62138300,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr18',
            'start': 59805532,
            'end': 59805533,
            'base': 'T'
        }
    },
    'chr19:55154216':{
        '38_coordinates': {
            'chrom': 'chr19',
            'start': 55154216,
            'end': 55154217,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr19',
            'start': 55665584,
            'end': 55665585,
            'base': 'A'
        }
    },
    'chrX:77682471':{
        '38_coordinates': {
            'chrom': 'chrX',
            'start': 77682471,
            'end': 77682472,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chrX',
            'start': 76937963,
            'end': 76937964,
            'base': 'G'
        }
    },
    'chrX:78043359': {
        '38_coordinates': {
            'chrom': 'chrX',
            'start': 78043359,
            'end': 78043360,
            'base': 'A'
        },
        '37_coordinates': {
            'chrom': 'chrX',
            'start': 77298857,
            'end': 77298858,
            'base': 'G'
        }
    },
}

ALZ_MISMATCH_SITES = {
    'chr2:21012603':{
        '38_coordinates': {
            'chrom': 'chr2',
            'start': 21012603,
            'end': 21012604,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr2',
            'start': 21235475,
            'end': 21235476,
            'base': 'T'
        }
    },
    'chr2:73448097':{
        '38_coordinates': {
            'chrom': 'chr2',
            'start': 73448097,
            'end': 73448101,
            'base': 'TCTC'
        },
        '37_coordinates': {
            'chrom': 'chr2',
            'start': 73675227,
            'end': 73675228,
            'base': 'T'
        }
    },
    'chr2:178756750':{
        '38_coordinates': {
            'chrom': 'chr2',
            'start': 178756750,
            'end': 178756751,
            'base': 'T'
        },
        '37_coordinates': {
            'chrom': 'chr2',
            'start': 179621477,
            'end': 179621478,
            'base': 'C'
        }
    },
    'chr6:7563750':{
        '38_coordinates': {
            'chrom': 'chr6',
            'start': 7563750,
            'end': 7563751,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr6',
            'start': 7563983,
            'end': 7563984,
            'base': 'T'
        }
    },
    'chr6:112187567':{
        '38_coordinates': {
            'chrom': 'chr6',
            'start': 112187567,
            'end': 112187568,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr6',
            'start': 112508769,
            'end': 112508770,
            'base': 'T'
        }
    },
    'chr6:112187568':{
        '38_coordinates': {
            'chrom': 'chr6',
            'start': 112187568,
            'end': 112187569,
            'base': 'T'
        },
        '37_coordinates': {
            'chrom': 'chr6',
            'start': 112508770,
            'end': 112508771,
            'base': 'G'
        }
    },
    'chr8:143213308':{
        '38_coordinates': {
            'chrom': 'chr8',
            'start': 143213308,
            'end': 143213309,
            'base': 'T'
        },
        '37_coordinates': {
            'chrom': 'chr8',
            'start': 144295183,
            'end': 144295184,
            'base': 'G'
        }
    },
    'chr10:110812700':{
        '38_coordinates': {
            'chrom': 'chr10',
            'start': 110812700,
            'end': 110812701,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr10',
            'start': 112572458,
            'end': 112572459,
            'base': 'G'
        }
    },
    'chr10:110831052':{
        '38_coordinates': {
            'chrom': 'chr10',
            'start': 110831052,
            'end': 110831053,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr10',
            'start': 112590810,
            'end': 112590811,
            'base': 'G'
        }
    },
    'chr11:121532520':{
        '38_coordinates': {
            'chrom': 'chr11',
            'start': 121532520,
            'end': 121532521,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr11',
            'start': 121403229,
            'end': 121403230,
            'base': 'T'
        }
    },
    'chr12:2682039':{
        '38_coordinates': {
            'chrom': 'chr12',
            'start': 2682039,
            'end': 2682040,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr12',
            'start': 2791205,
            'end': 2791206,
            'base': 'A'
        }
    },
    'chr12:21864488':{
        '38_coordinates': {
            'chrom': 'chr12',
            'start': 21864488,
            'end': 21864489,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr12',
            'start': 22017422,
            'end': 22017423,
            'base': 'A'
        }
    },
    'chr12:21910181':{
        '38_coordinates': {
            'chrom': 'chr12',
            'start': 21910181,
            'end': 21910182,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr12',
            'start': 22063115,
            'end': 22063116,
            'base': 'A'
        }
    },
    'chr13:110187203':{
        '38_coordinates': {
            'chrom': 'chr13',
            'start': 110187203,
            'end': 110187204,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr13',
            'start': 110839550,
            'end': 110839551,
            'base': 'T'
        }
    },
    'chr15:48515440':{
        '38_coordinates': {
            'chrom': 'chr15',
            'start': 48515440,
            'end': 48515441,
            'base': 'T'
        },
        '37_coordinates': {
            'chrom': 'chr15',
            'start': 48807637,
            'end': 48807638,
            'base': 'C'
        }
    },
    'chr17:50188134':{
        '38_coordinates': {
            'chrom': 'chr17',
            'start': 50188134,
            'end': 50188135,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr17',
            'start': 48265495,
            'end': 48265496,
            'base': 'T'
        }
    },
    'chr17:80325138':{
        '38_coordinates': {
            'chrom': 'chr17',
            'start': 80325138,
            'end': 80325139,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr17',
            'start': 78298938,
            'end': 78298939,
            'base': 'A'
        }
    },
    'chr19:40617187':{
        '38_coordinates': {
            'chrom': 'chr19',
            'start': 40617187,
            'end': 40617189,
            'base': 'AG'
        },
        '37_coordinates': {
            'chrom': 'chr19',
            'start': 41123093,
            'end': 41123094,
            'base': 'A'
        }
    },
    'chr19:55154216':{
        '38_coordinates': {
            'chrom': 'chr19',
            'start': 55154216,
            'end': 55154217,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr19',
            'start': 55665584,
            'end': 55665585,
            'base': 'A'
        }
    },
    'chr21:14144655':{
        '38_coordinates': {
            'chrom': 'chr21',
            'start': 14144655,
            'end': 14144656,
            'base': 'A'
        },
        '37_coordinates': {
            'chrom': 'chr21',
            'start': 15516976,
            'end': 15516977,
            'base': 'T'
        }
    },
}


CARDIOMET_MISMATCH_SITES = {
    'chr2:21012603':{
        '38_coordinates': {
            'chrom': 'chr2',
            'start': 21012603,
            'end': 21012604,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr2',
            'start': 21235475,
            'end': 21235476,
            'base': 'T'
        }
    },
    'chr2:73448097':{
        '38_coordinates': {
            'chrom': 'chr2',
            'start': 73448097,
            'end': 73448101,
            'base': 'TCTC'
        },
        '37_coordinates': {
            'chrom': 'chr2',
            'start': 73675227,
            'end': 73675228,
            'base': 'T'
        }
    },
    'ch2:178756750':{
        '38_coordinates': {
            'chrom': 'chr2',
            'start': 178756750,
            'end': 178756751,
            'base': 'T'
        },
        '37_coordinates': {
            'chrom': 'chr2',
            'start': 179621477,
            'end': 179621478,
            'base': 'C'
        }
    },
    'chr6:7563750':{
        '38_coordinates': {
            'chrom': 'chr6',
            'start': 7563750,
            'end': 7563751,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr6',
            'start': 7563983,
            'end': 7563984,
            'base': 'T'
        }
    },
    'chr6:112187567':{
        '38_coordinates': {
            'chrom': 'chr6',
            'start': 112187567,
            'end': 112187568,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr6',
            'start': 112508769,
            'end': 112508770,
            'base': 'T'
        }
    },
    'chr6: 112187568':{
        '38_coordinates': {
            'chrom': 'chr6',
            'start': 112187568,
            'end': 112187569,
            'base': 'T'
        },
        '37_coordinates': {
            'chrom': 'chr6',
            'start': 112508770,
            'end': 112508771,
            'base': 'G'
        }
    },
    'chr8:143213308':{
        '38_coordinates': {
            'chrom': 'chr8',
            'start': 143213308,
            'end': 143213309,
            'base': 'T'
        },
        '37_coordinates': {
            'chrom': 'chr8',
            'start': 144295183,
            'end': 144295184,
            'base': 'G'
        }
    },
    'chr10:110812700':{
        '38_coordinates': {
            'chrom': 'chr10',
            'start': 110812700,
            'end': 110812701,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr10',
            'start': 112572458,
            'end': 112572459,
            'base': 'G'
        }
    },
    'chr10:110831052':{
        '38_coordinates': {
            'chrom': 'chr10',
            'start': 110831052,
            'end': 110831053,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr10',
            'start': 112590810,
            'end': 112590811,
            'base': 'G'
        }
    },
    'chr12:2682039':{
        '38_coordinates': {
            'chrom': 'chr12',
            'start': 2682039,
            'end': 2682040,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr12',
            'start': 2791205,
            'end': 2791206,
            'base': 'A'
        }
    },
    'chr12:21910181':{
        '38_coordinates': {
            'chrom': 'chr12',
            'start': 21910181,
            'end': 21910182,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr12',
            'start': 22063115,
            'end': 22063116,
            'base': 'A'
        }
    },
    'chr12:21864488':{
        '38_coordinates': {
            'chrom': 'chr12',
            'start': 21864488,
            'end': 21864489,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr12',
            'start': 22017422,
            'end': 22017423,
            'base': 'A'
        }
    },
    'chr13:110187203':{
        '38_coordinates': {
            'chrom': 'chr13',
            'start': 110187203,
            'end': 110187204,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr13',
            'start': 110839550,
            'end': 110839551,
            'base': 'T'
        }
    },
    'chr15:48515440':{
        '38_coordinates': {
            'chrom': 'chr15',
            'start': 48515440,
            'end': 48515441,
            'base': 'T'
        },
        '37_coordinates': {
            'chrom': 'chr15',
            'start': 48807637,
            'end': 48807638,
            'base': 'C'
        }
    },
    'chr17:50188134':{
        '38_coordinates': {
            'chrom': 'chr17',
            'start': 50188134,
            'end': 50188135,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr17',
            'start': 48265495,
            'end': 48265496,
            'base': 'T'
        }
    },
    'chr17:80325138':{
        '38_coordinates': {
            'chrom': 'chr17',
            'start': 80325138,
            'end': 80325139,
            'base': 'G'
        },
        '37_coordinates': {
            'chrom': 'chr17',
            'start': 78298938,
            'end': 78298939,
            'base': 'A'
        }
    },
    'chr19:40617187':{
        '38_coordinates': {
            'chrom': 'chr19',
            'start': 40617187,
            'end': 40617189,
            'base': 'AG'
        },
        '37_coordinates': {
            'chrom': 'chr19',
            'start': 41123093,
            'end': 41123094,
            'base': 'A'
        }
    },
    'chr19:55154216':{
        '38_coordinates': {
            'chrom': 'chr19',
            'start': 55154216,
            'end': 55154217,
            'base': 'C'
        },
        '37_coordinates': {
            'chrom': 'chr19',
            'start': 55665584,
            'end': 55665585,
            'base': 'A'
        }
    },
    'chr21:14144655':{
        '38_coordinates': {
            'chrom': 'chr21',
            'start': 14144655,
            'end': 14144656,
            'base': 'A'
        },
        '37_coordinates': {
            'chrom': 'chr21',
            'start': 15516976,
            'end': 15516977,
            'base': 'T'
        }
    },
}


def find_overlapping_mismatch_site(record,MISMATCH_SITES):
    for key, site in MISMATCH_SITES.items():
        if record.CHROM == site['38_coordinates']['chrom']:
            if record.POS <= site['38_coordinates']['start'] <= record.end:
                return key


def update_grch38_ref_to_grch37_for_record_if_needed(record, MISMATCH_SITES, mismatched_site_key=None):
    """
    If record overlaps mismatched sites, update ref, alt and genotype accordingly
    Some assumptions:
    - records have only one sample (true for color data, one file per sample)
    - will not try to fix record, if genotype is malformed or missing
    - will not try to fix record, if more than 2 alleles
    """
    if not mismatched_site_key:
        mismatched_site_key = find_overlapping_mismatch_site(record, MISMATCH_SITES)
    if not mismatched_site_key:
        return record
    mismatched_site = MISMATCH_SITES[mismatched_site_key]
    logger = logging.getLogger(__name__)
    # get genotype
    # code only for one sample per vcf file. This is how color files are
    assert (len(record.samples) == 1)
    sample_call = record.samples[0]
    if 'GT' not in record.FORMAT:
        # missing genotype issue not handled
        raise ValueError(f'GT not in record {record}')
    gt_indices = sample_call.gt_alleles
    expected_gts = set(str(i) for i in range(len(record.ALT) + 1))
    if not set(gt_indices).issubset(expected_gts):
        # unknown genotype issue not handled
        raise ValueError(f'Unknown genotypes for record {record}, genotype {gt_indices}')

    ref = list(record.REF)
    pos = mismatched_site['38_coordinates']['start'] - record.POS

    if mismatched_site['38_coordinates']['base'] == ref[pos]:
        ref[pos] = mismatched_site['37_coordinates']['base']
    updated_ref = ''.join(ref)
    orig_allele_seqs = []
    # alleles from original record
    for allele in sample_call.gt_alleles:
        # only add reference, if it is in genotype, het variant
        if allele == '0':
            orig_allele_seqs.append(record.REF)
        else:
            orig_allele_seqs.append(record.ALT[int(allele) - 1])
    updated_alts = []
    updated_gt = []
    for allele_seq in orig_allele_seqs:
        # original allele matches updated ref, genotype needs '0' added
        if allele_seq == updated_ref:
            updated_gt.append(0)
        else:
            if allele_seq not in updated_alts:
                updated_alts.append(allele_seq)
            updated_gt.append(updated_alts.index(allele_seq) + 1)
    updated_gt = '/'.join(map(str, sorted(updated_gt)))

    if len(updated_alts) > 2:
        # case we will not try to handle, write original record as is
        raise ValueError(f'Updated record has more than 2 alts {updated_alts}, original record {record}')

    if len(updated_alts) == 0:
        # return None if the variant was not a variant in 37
        record = None
    else:
        # add genotype data to record. Loses all other fields
        record.FORMAT = 'GT'
        record.samples[0].data = calldata_spec(updated_gt)

        record.REF = updated_ref
        record.ALT = []
        record.add_info('PREPROCESSED')
        for alt in updated_alts:
            record.ALT.append(vcf.model._Substitution(alt))

    return record


def _create_fake_record():
    """
    If the input vcf file is empty, use a fake record to use as template to create mismatch records
    """
    with StringIO('\n'.join(EXAMPLE_VCF)) as fh:
        records = list(vcf.Reader(fh))
        return records[0]


def convert_grch38_ref_mismatch_sites_to_grch37(input_vcf_file, output_vcf_basename,project):
    """
    For ACMG59 reportable range there are 4 sites that have
    reference mismatch between GRCh37 and GRCh38
    All ref and alts in variants overlapping these sites
    will need to be updated to 37 reference
    output file will contain variants overlapping mismatch sites
    and all other variants with original record
    """
    logger = logging.getLogger(__name__)
    output_vcf_file = f'{output_vcf_basename}.vcf'
    reader = vcf.Reader(filename=input_vcf_file)
    records = list(reader)
    if len(records) == 0:
        record = _create_fake_record()
        logger.warning(f'VCF file {input_vcf_file} is empty. We will be creating the mismatch sites using a template')

    mismatched_site_overlap = {}
    if project == 'aou':
            MISMATCH_SITES = AOU_MISMATCH_SITES
            print('Using: AOU_MISMATCH_SITES')
    elif project == 'epilepsy':
            MISMATCH_SITES = EPILEPSY_MISMATCH_SITES
            print('Using: EPILEPSY_MISMATCH_SITES')
    elif project == 'alzheimers':
            MISMATCH_SITES = ALZ_MISMATCH_SITES
            print('Using: ALZ_MISMATCH_SITES')
    elif project == 'cardiomet':
            MISMATCH_SITES = CARDIOMET_MISMATCH_SITES
            print('Using: CARDIOMET_MISMATCH_SITES')

    for record in records:
        mismatched_site_key = find_overlapping_mismatch_site(record,MISMATCH_SITES)
        if mismatched_site_key:
            mismatched_site_overlap[mismatched_site_key] = True
            try:
                update_grch38_ref_to_grch37_for_record_if_needed(record, MISMATCH_SITES, mismatched_site_key)
            except ValueError as e:
                logger.info(f'Record {record.CHROM}:{record.POS} with mismatch site {mismatched_site_key} encountered error {e}')

    reader.infos['PREPROCESSED'] = VcfInfo(
        'PREPROCESSED',
        0,
        'Flag',
        'The record was pre-processed. Added when a record needed to be changed for liftover',
        '',
        '',
    )
    # if there are no overlapping variants in mismatched sites,
    # create a homozygous variant matching 37 as ref and 38 as alt
    for key, site in MISMATCH_SITES.items():
        if key not in mismatched_site_overlap.keys():
            # TODO: separate out creation of a record
            mismatch_record = copy(record)
            mismatch_record.ID = '.'
            mismatch_record.QUAL = 100
            mismatch_record.FILTER = []
            mismatch_record.FORMAT = 'GT'
            mismatch_record.samples = []
            # copy the objects within a record.
            # Without doing an explicit copy it will just be a
            # pointer to the original record
            for sample in record.samples:
                mismatch_record.samples.append(copy(sample))
            mismatch_record.samples[0].data = calldata_spec('1/1')
            mismatch_record.INFO = {}
            mismatch_record.add_info('PREPROCESSED')
            mismatch_record.CHROM = site['38_coordinates']['chrom']
            mismatch_record.POS = site['38_coordinates']['start']
            mismatch_record.REF = site['37_coordinates']['base']
            mismatch_record.ALT = [
                vcf.model._Substitution(site['38_coordinates']['base'])
            ]
            records.append(mismatch_record)

    contig_order = {c: i for i, c in enumerate(reader.contigs)}

    def sort_key(record):
        """
        Sorts records by (CHROM,POS,REF).
        If contigs are specified in the VCF file and record CHROM matches a contig,
        contig order is maintained.
        Any unmatched CHROMs will throw an error
        """
        if record.CHROM not in contig_order:
            raise ValueError(
                f'Unexpected chrom {record.CHROM} found. Expected one of {contig_order.keys()}'
            )
        return (contig_order[record.CHROM], record.POS, record.REF)

    records.sort(key=sort_key)

    with open(output_vcf_file, 'w') as out_fp:
        writer = vcf.Writer(out_fp, reader, lineterminator='\n')
        for record in records:
            writer.write_record(record)



if __name__ == '__main__':
    function = getattr(sys.modules[__name__], sys.argv[1])
    function(sys.argv[2], sys.argv[3], sys.argv[4])