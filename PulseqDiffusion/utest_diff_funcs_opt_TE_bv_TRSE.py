import unittest
import diff_funcs as difunc
import pickle
import sys

from pypulseq.Sequence.sequence import Sequence

class TestOptTEbvTRSE(unittest.TestCase):
    def test_optTEbvTRSE1(self):
        with open('utest_sample_data/opt_TE_bv_TRSE1.pkl', 'rb') as f:
            bvalue_Dict, grads_times_Dict, seq_sys_Dict, result_tmp = pickle.load(f)

        result = difunc.opt_TE_bv_TRSE(bvalue_Dict, grads_times_Dict, seq_sys_Dict)
        self.assertEqual(result[0:2], result_tmp[0:2])
        self.assertTrue(self, all(result[2] == result_tmp[2]))
        self.assertEqual(result[3:], result_tmp[3:])
        print('------End test 1------')

    def test_optTEbvTRSE2(self):
        with open('utest_sample_data/opt_TE_bv_TRSE2.pkl', 'rb') as f:
            bvalue_Dict, grads_times_Dict, seq_sys_Dict, result_tmp = pickle.load(f)

        result = difunc.opt_TE_bv_TRSE(bvalue_Dict, grads_times_Dict, seq_sys_Dict)
        self.assertEqual(result[0:2], result_tmp[0:2])
        self.assertTrue(self, all(result[2] == result_tmp[2]))
        self.assertEqual(result[3:], result_tmp[3:])
        print('------End test 2------')

    def test_optTEbvTRSE3(self):
        with open('utest_sample_data/opt_TE_bv_TRSE3.pkl', 'rb') as f:
            bvalue_Dict, grads_times_Dict, seq_sys_Dict, result_tmp = pickle.load(f)

        result = difunc.opt_TE_bv_TRSE(bvalue_Dict, grads_times_Dict, seq_sys_Dict)
        self.assertEqual(result[0:2], result_tmp[0:2])
        self.assertTrue(self, all(result[2] == result_tmp[2]))
        self.assertEqual(result[3:], result_tmp[3:])
        print('------End test 3------')

if __name__ == '__main__':
    unittest.main()