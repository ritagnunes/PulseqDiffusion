import unittest
import diff_funcs as difunc

class TestCalbval(unittest.TestCase):
    def test_calbval1(self):
        max_grad = 1362432.0 # Hz
        n_gdiff_delta = 1716
        i_raster_time = 100000
        n_gdiff_Delta = 2892
        n_gdiff_rt = 22
        result = difunc.calc_bval(max_grad, n_gdiff_delta / i_raster_time, n_gdiff_Delta / i_raster_time, n_gdiff_rt / i_raster_time)
        self.assertEqual(result, 500614220.8374321)

    def test_calbval2(self):
        max_grad = 1362432.0  # Hz
        n_gdiff_delta = 2257
        i_raster_time = 100000
        n_gdiff_Delta = 3433
        n_gdiff_rt = 22
        result = difunc.calc_bval(max_grad, n_gdiff_delta / i_raster_time, n_gdiff_Delta / i_raster_time, n_gdiff_rt / i_raster_time)
        self.assertEqual(result,1000666960.8932015)

    def test_calbval3(self):
        max_grad = 1362432.0  # Hz
        n_gdiff_delta = 1000
        i_raster_time = 100000
        n_gdiff_Delta = 2000
        n_gdiff_rt = 10
        result = difunc.calc_bval(max_grad, n_gdiff_delta / i_raster_time, n_gdiff_Delta / i_raster_time, n_gdiff_rt / i_raster_time)
        self.assertEqual(result, 122133224.45259748)


if __name__ == '__main__':
    unittest.main()
