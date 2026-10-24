import csv
import io

COLUMNS = 'proband_id chrom_1 pos_1 ref_1 alt_1 chrom_2 pos_2 ref_2 alt_2 epcr finding_type notes'.split()

def test_synthetic_schema_and_order():
    synthetic = ','.join(COLUMNS) + '\nSYN,chr1,2,C,T,,,,,0.1,primary,synthetic\nSYN,chr1,3,G,A,,,,,0.05,primary,synthetic\n'
    reader = csv.DictReader(io.StringIO(synthetic)); rows = list(reader)
    assert reader.fieldnames == COLUMNS and len(rows) == 2
    assert [float(x['epcr']) for x in rows] == sorted([float(x['epcr']) for x in rows], reverse=True)

def test_synthetic_ref_and_ties():
    reference = 'ACGT'
    assert reference[1] == 'C'
    assert sorted([(0.1, 'B'), (0.1, 'A')], key=lambda x:(-x[0],x[1])) == [(0.1,'A'),(0.1,'B')]
