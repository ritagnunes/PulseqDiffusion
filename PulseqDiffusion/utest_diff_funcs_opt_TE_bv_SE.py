import unittest
import diff_funcs as difunc
import pickle
import sys

from pypulseq.Sequence.sequence import Sequence

class TestOptTEbvSE(unittest.TestCase):
    def test_optTEbvSE1(self):
        with open('utest_sample_data/opt_TE_bv_SE1.pkl', 'rb') as f:
            bvalue_Dict, grads_times_Dict, seq_sys_Dict, result_tmp = pickle.load(f)

        result = difunc.opt_TE_bv_SE(bvalue_Dict, grads_times_Dict, seq_sys_Dict)
        self.assertEqual(result, result_tmp)
        print('------End test 1------')

    def test_optTEbvSE2(self):
        with open('utest_sample_data/opt_TE_bv_SE2.pkl', 'rb') as f:
            bvalue_Dict, grads_times_Dict, seq_sys_Dict, result_tmp = pickle.load(f)

        result = difunc.opt_TE_bv_SE(bvalue_Dict, grads_times_Dict, seq_sys_Dict)
        self.assertEqual(result, result_tmp)
        print('------End test 2------')

    def test_optTEbvSE3(self):
        with open('utest_sample_data/opt_TE_bv_SE3.pkl', 'rb') as f:
            bvalue_Dict, grads_times_Dict, seq_sys_Dict, result_tmp = pickle.load(f)

        result = difunc.opt_TE_bv_SE(bvalue_Dict, grads_times_Dict, seq_sys_Dict)
        self.assertEqual(result, result_tmp)
        print('------End test 3------')

if __name__ == '__main__':
    unittest.main()