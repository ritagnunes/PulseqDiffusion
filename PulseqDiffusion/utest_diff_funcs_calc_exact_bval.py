import unittest
import diff_funcs as difunc
import pickle
import sys

from pypulseq.Sequence.sequence import Sequence

class TestCalexactbval(unittest.TestCase):
    def test_calexacetbval1(self):
        with open('utest_sample_data/calc_exact_bval1.pkl', 'rb') as f:
            waveform, INV, short_f, seq = pickle.load(f)

        result = difunc.calc_exact_bval(waveform, INV, short_f, seq)
        self.assertEqual(result, 327.19201047729877)

    def test_calexacetbval2(self):
        with open('utest_sample_data/calc_exact_bval2.pkl', 'rb') as f:
            waveform, INV, short_f, seq = pickle.load(f)

        result = difunc.calc_exact_bval(waveform, INV, short_f, seq)
        self.assertEqual(result, 498.775665378322)

    def test_calexacetbval3(self):
        with open('utest_sample_data/calc_exact_bval3.pkl', 'rb') as f:
            waveform, INV, short_f, seq = pickle.load(f)

        result = difunc.calc_exact_bval(waveform, INV, short_f, seq)
        self.assertEqual(result, 749.6658744624788)

if __name__ == '__main__':
    unittest.main()