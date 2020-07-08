import pickle
import sys
import unittest
from pathlib import Path

from pypulseq import opts
from pypulseq.Sequence import sequence

import PulseqDiffusion.diff_funcs as difunc

sys.modules['sequence'] = sequence
sys.modules['opts'] = opts


class TestCalexactbval(unittest.TestCase):
    def test_calexacetbval1(self):
        path = Path(__file__).parent.parent / 'utest_sample_data' / 'calc_exact_bval1.pkl'
        with open(path, 'rb') as f:
            waveform, INV, short_f, seq = pickle.load(f)

        result = difunc.calc_exact_bval(waveform, INV, short_f, seq)
        self.assertEqual(result, 327.19201047729877)

    def test_calexacetbval2(self):
        path = Path(__file__).parent.parent / 'utest_sample_data' / 'calc_exact_bval2.pkl'
        with open(path, 'rb') as f:
            waveform, INV, short_f, seq = pickle.load(f)

        result = difunc.calc_exact_bval(waveform, INV, short_f, seq)
        self.assertEqual(result, 498.775665378322)

    def test_calexacetbval3(self):
        path = Path(__file__).parent.parent / 'utest_sample_data' / 'calc_exact_bval3.pkl'
        with open(path, 'rb') as f:
            waveform, INV, short_f, seq = pickle.load(f)

        result = difunc.calc_exact_bval(waveform, INV, short_f, seq)
        self.assertEqual(result, 749.6658744624788)


if __name__ == '__main__':
    unittest.main()
