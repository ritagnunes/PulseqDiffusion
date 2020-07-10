import pickle
import sys
import unittest
from pathlib import Path

from pypulseq import opts
from pypulseq.Sequence import sequence

if __name__ == '__main__':
    path = Path(__file__).absolute().parent.parent.parent
    sys.path.insert(0, str(path))

import PulseqDiffusion.diff_funcs as difunc

sys.modules['sequence'] = sequence
sys.modules['opts'] = opts


class TestOptTEbvTRSE(unittest.TestCase):
    def test_optTEbvTRSE1(self):
        path = Path(__file__).parent / 'utest_sample_data' / 'opt_TE_bv_TRSE1.pkl'
        with open(path, 'rb') as f:
            bvalue_Dict, grads_times_Dict, seq_sys_Dict, result_tmp = pickle.load(f)

        result = difunc.opt_TE_bv_TRSE(bvalue_Dict, grads_times_Dict, seq_sys_Dict)
        self.assertAlmostEqual(result[0], result_tmp[0])
        self.assertAlmostEqual(result[1], result_tmp[1])
        self.assertTrue(self, all(result[2] == result_tmp[2]))
        self.assertAlmostEqual(result[3:], result_tmp[3:])
        print('------End test 1------')

    def test_optTEbvTRSE2(self):
        path = Path(__file__).parent / 'utest_sample_data' / 'opt_TE_bv_TRSE2.pkl'
        with open(path, 'rb') as f:
            bvalue_Dict, grads_times_Dict, seq_sys_Dict, result_tmp = pickle.load(f)

        result = difunc.opt_TE_bv_TRSE(bvalue_Dict, grads_times_Dict, seq_sys_Dict)
        self.assertAlmostEqual(result[0], result_tmp[0])
        self.assertAlmostEqual(result[1], result_tmp[1])
        self.assertTrue(self, all(result[2] == result_tmp[2]))
        self.assertAlmostEqual(result[3:], result_tmp[3:])
        print('------End test 2------')

    def test_optTEbvTRSE3(self):
        path = Path(__file__).parent / 'utest_sample_data' / 'opt_TE_bv_TRSE3.pkl'
        with open(path, 'rb') as f:
            bvalue_Dict, grads_times_Dict, seq_sys_Dict, result_tmp = pickle.load(f)

        result = difunc.opt_TE_bv_TRSE(bvalue_Dict, grads_times_Dict, seq_sys_Dict)
        self.assertAlmostEqual(result[0], result_tmp[0])
        self.assertAlmostEqual(result[1], result_tmp[1])
        self.assertTrue(self, all(result[2] == result_tmp[2]))
        self.assertAlmostEqual(result[3:], result_tmp[3:])
        print('------End test 3------')


if __name__ == '__main__':
    unittest.main()
