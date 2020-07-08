import unittest

from PulseqDiffusion.test.utest_diff_funcs_calc_bval import TestCalbval
from PulseqDiffusion.test.utest_diff_funcs_calc_exact_bval import TestCalexactbval
from PulseqDiffusion.test.utest_diff_funcs_getdirs import TestGetDirs
from PulseqDiffusion.test.utest_diff_funcs_opt_TE_bv_SE import TestOptTEbvSE
from PulseqDiffusion.test.utest_diff_funcs_opt_TE_bv_TRSE import TestOptTEbvTRSE


def create_suite():
    test_suite = unittest.TestSuite()
    test_suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestCalbval))
    test_suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestCalexactbval))
    test_suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestGetDirs))
    test_suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestOptTEbvSE))
    test_suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestOptTEbvTRSE))
    return test_suite


if __name__ == '__main__':
    suite = create_suite()
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
    report = f'{result.testsRun} tests run\n'
    report += f'All tests successful: {result.wasSuccessful()}\n'
    report += f'Errors: {result.errors}'
    with open('Unit test result.txt', 'w') as f:
        f.write(report)
