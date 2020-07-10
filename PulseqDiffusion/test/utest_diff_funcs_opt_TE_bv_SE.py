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


class TestOptTEbvSE(unittest.TestCase):
    def test_optTEbvSE1(self):
        path = Path(__file__).parent / 'utest_sample_data' / 'opt_TE_bv_SE1.pkl'
        with open(path, 'rb') as f:
            bvalue_Dict, grads_times_Dict, seq_sys_Dict, result_tmp = pickle.load(f)

        result = difunc.opt_TE_bv_SE(bvalue_Dict, grads_times_Dict, seq_sys_Dict)
        self.assertEqual(result, result_tmp)
        print('------End test 1------')

    def test_optTEbvSE2(self):
        path = Path(__file__).parent / 'utest_sample_data' / 'opt_TE_bv_SE2.pkl'
        with open(path, 'rb') as f:
            bvalue_Dict, grads_times_Dict, seq_sys_Dict, result_tmp = pickle.load(f)

        result = difunc.opt_TE_bv_SE(bvalue_Dict, grads_times_Dict, seq_sys_Dict)
        self.assertEqual(result, result_tmp)
        print('------End test 2------')

    def test_optTEbvSE3(self):
        path = Path(__file__).parent/ 'utest_sample_data' / 'opt_TE_bv_SE3.pkl'
        with open(path, 'rb') as f:
            bvalue_Dict, grads_times_Dict, seq_sys_Dict, result_tmp = pickle.load(f)

        result = difunc.opt_TE_bv_SE(bvalue_Dict, grads_times_Dict, seq_sys_Dict)
        self.assertEqual(result, result_tmp)
        print('------End test 3------')


if __name__ == '__main__':
    unittest.main()
