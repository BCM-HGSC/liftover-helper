import collections
import os
import tempfile
import unittest

import vcf

from scripts import liftover

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


class RecordTest(unittest.TestCase):
    def test_convert_hg19_vcf_to_grch37_vcf(self):
        temp_dir = tempfile.gettempdir()
        output_file = os.path.join(temp_dir, 'output_grch37.vcf')
        liftover.convert_hg19_vcf_to_grch37_vcf('tests/hg19.vcf', output_file)

        with open(output_file) as fh:
            records = vcf.Reader(fh)
            expected_contigs = collections.OrderedDict([('1', liftover.contig_spec('1', 249250621)),
                                                        ('2', liftover.contig_spec('2', 243199373)),
                                                        ('3', liftover.contig_spec('3', 198022430))])
            self.assertEqual(records.contigs, expected_contigs)
            record = next(records)
            self.assertEqual(record.CHROM, '1')
            self.assertEqual(record.POS, 97915604)
            record = next(records)
            self.assertEqual(record.CHROM, '2')
            self.assertEqual(record.POS, 97915605)


    def test_record_overlaps_mismatch(self):
        records = list(vcf.Reader(filename='tests/grch38.vcf'))
        expected_result = 'chr2:21012603'
        for i in range(0, 5):
            self.assertEqual(liftover.find_overlapping_mismatch_site(records[i], AOU_MISMATCH_SITES), expected_result)
        for i in range(5, 7):
            self.assertFalse(liftover.find_overlapping_mismatch_site(records[i], AOU_MISMATCH_SITES))


    def test_update_record(self):
        records = list(vcf.Reader(filename='tests/grch38_2.vcf'))
        expected_results = [
            ['ATATG', 'ACATG,A', '1/2'],
            ['ATATG', 'A', '1/1'],
            ['AT', 'AC,A', '1/2'],
            ['AT', 'A', '1/1'],
            ['T', 'C,A', '1/2'],
            ['T', 'A', '1/1'],
            ['T', 'C', '0/1'],
            None,
            ['C', 'T', './.', True],
            ['C', 'T', '1/2', True],
            ['TATG', 'CATG,C', '1/2'],
            ['TATG', 'C', '1/1'],
            ['T', 'C,CAAT', '1/2'],
            ['T', 'CAAT', '0/1'],
            ['ATG', 'A', '0/1'],
            ['C', 'T', '0/1'],
        ]
        for i in range(len(records)):
            expected_record = expected_results[i]
            if expected_record and len(expected_record) == 4:
                self.assertRaises(ValueError, liftover.update_grch38_ref_to_grch37_for_record_if_needed, records[i], AOU_MISMATCH_SITES)
                continue

            observed_record = liftover.update_grch38_ref_to_grch37_for_record_if_needed(records[i], AOU_MISMATCH_SITES)
            if expected_record is None:
                self.assertIsNone(observed_record)
            else:
                self.assertEqual(observed_record.REF, expected_record[0])
                observed_alts = ','.join(map(str, observed_record.ALT))
                self.assertEqual(observed_alts, expected_record[1])
                self.assertEqual(observed_record.samples[0].data.GT, expected_record[2])

        # TODO: tests for is_anchor_base

if __name__ == '__main__':
    unittest.main()
