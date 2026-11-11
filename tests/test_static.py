"""Run with python -m unittest discover -s tests; requires no Hail install."""
import ast
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts'

class TestSource(unittest.TestCase):
    def test_python_syntax(self):
        for p in SCRIPTS.glob('*.py'):
            with self.subTest(p=p.name):
                ast.parse(p.read_text())

    def test_student_filtering_block_identical(self):
        orig = (SCRIPTS / 'student_original.py').read_text()
        newer = (SCRIPTS / '02_mendelian_screen.py').read_text()
        a = '    # Recompute original *trio-cohort, pre-QC* AC'
        b = '    # Counts refer to the final TSV/VCF'
        self.assertEqual(orig[orig.index(a):orig.index(b)], newer[newer.index(a):newer.index(b)])

    def test_screen_reads_prepared_matrix(self):
        content = (SCRIPTS / '02_mendelian_screen.py').read_text()
        self.assertIn('hl.read_matrix_table(args.matrix)', content)
        self.assertNotIn('hl.import_vcf(', content)

    def test_prepare_writes_durable_matrix(self):
        content = (SCRIPTS / '01_prepare_matrix.py').read_text()
        self.assertIn('mt.write(a.matrix, overwrite=False)', content)
        self.assertIn('validate_existing(a)', content)

if __name__ == '__main__':
    unittest.main()
