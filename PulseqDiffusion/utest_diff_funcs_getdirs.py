import unittest
import numpy as np
import diff_funcs as difunc

class TestGetDirs(unittest.TestCase):
    def test_getdirs1(self):
        ndirs=3
        gutest, nb0s = difunc.get_dirs(ndirs)
        gpred = np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
        self.assertEqual(nb0s, 1)
        np.testing.assert_allclose(gpred, gutest)

    def test_getdirs2(self):
        ndirs=6
        gutest, nb0s = difunc.get_dirs(ndirs)
        gpred = np.array([[-0.283341,  -0.893706,  -0.347862],
                     [-0.434044,   0.799575,  -0.415074],
                     [0.961905,    0.095774,   0.256058],
                     [-0.663896,   0.491506,   0.563616],
                     [-0.570757,  -0.554998,   0.605156],
                     [-0.198848,  -0.056534,  -0.978399]])
        self.assertEqual(nb0s, 1)
        np.testing.assert_allclose(gpred, gutest)

    def test_getdirs3(self):
        ndirs=12
        gutest, nb0s = difunc.get_dirs(ndirs)
        gpred = np.array([[ 0.648514,   0.375307,   0.662249],
                     [-0.560493,   0.824711,   0.075496],
                     [-0.591977,  -0.304283,   0.746307],
                     [-0.084472,  -0.976168,   0.199902],
                     [-0.149626,  -0.006494,  -0.988721],
                     [-0.988211,  -0.056904,  -0.142130],
                     [0.864451,   -0.274379,  -0.421237],
                     [-0.173549,  -0.858586,  -0.482401],
                     [-0.039288,   0.644885,  -0.763269],
                     [0.729809,    0.673235,  -0.118882],
                     [0.698325,   -0.455759,   0.551929],
                     [-0.325340,   0.489774,   0.808873]])
        self.assertEqual(nb0s, 1)
        np.testing.assert_allclose(gpred, gutest)

    def test_getdirs4(self):
        ndirs=60
        gutest, nb0s = difunc.get_dirs(ndirs)
        gpred = np.array([[-0.811556,  0.245996, -0.529964],
                      [-0.576784, -0.313126,  0.754502],
                      [-0.167946, -0.899364, -0.403655],
                      [0.755699, -0.512113,     -0.408238],
                      [0.116846,   0.962654,    -0.244221],
                      [0.495465,   0.208081,    0.843337],
                      [0.901459,   0.385831,    -0.196230],
                      [-0.248754,   0.420519,   -0.872516],
                      [-0.047525,   0.444671,   0.894432],
                      [- 0.508593,   0.857494, -0.077699],
                      [0.693558,    0.614042,    0.376737],
                      [0.990394, -0.134781,    -0.030898],
                      [0.019140, -0.684235,     0.729010],
                      [0.385221, - 0.339346, - 0.858166],
                      [- 0.440289, - 0.853536,   0.278609],
                      [0.680515, - 0.559825,   0.472752],
                      [- 0.146029,   0.872237,   0.466774],
                      [0.317352,   0.195118, - 0.928018],
                      [- 0.796280, - 0.129004, - 0.591013],
                      [- 0.711299,   0.249255,   0.657211],
                      [- 0.838383, - 0.538321, - 0.085587],
                      [0.202544, - 0.966710,   0.156357],
                      [-0.296747, - 0.476761, - 0.827430],
                      [0.545225,   0.637023, - 0.544914],
                      [-0.887097,   0.451265,   0.097048],
                      [0.034752, - 0.124211,   0.991647],
                      [0.469222, - 0.766720, - 0.438145],
                      [-0.948457, - 0.088803,   0.304209],
                      [-0.354311,   0.664176, - 0.658281],
                      [-0.462117, - 0.833550, - 0.302724],
                      [0.949202, - 0.022353,  0.313871],
                      [0.791248,   0.065354, - 0.607994],
                      [-0.004026,   0.992213,   0.124488],
                      [- 0.357034,   0.107223, - 0.927917],
                      [0.414504,   0.685422,   0.598651],
                      [-0.331743, - 0.552720,   0.764492],
                      [-0.749058,   0.641807, - 0.164305],
                      [0.238666, - 0.655691, - 0.716315],
                      [0.619125,   0.784982,   0.022082],
                      [0.123966, - 0.872070,   0.473419],
                      [-0.185240,   0.122014,   0.975089],
                      [-0.980282, - 0.189151, - 0.057166],
                      [0.637873, - 0.084335,   0.765510],
                      [-0.668960,   0.723273,   0.171375],
                      [-0.775822, - 0.381543,   0.502519],
                      [-0.636044, - 0.425049, - 0.644035],
                      [0.229220,   0.809364, - 0.540729],
                      [-0.538340,   0.531550,   0.653946],
                      [0.906105, - 0.354129,   0.231445],
                      [-0.166743, - 0.191681, - 0.967189],
                      [0.324636, - 0.927784, - 0.183922],
                      [0.551291, - 0.398459,   0.733014],
                      [0.537753, - 0.032690, - 0.842468],
                      [-0.306182, - 0.951456, - 0.031381],
                      [0.875976,   0.329209,   0.352545],
                      [0.902989, - 0.218632, - 0.369879],
                      [-0.456427,   0.801551, - 0.386251],
                      [0.089001,   0.716134,   0.692265],
                      [-0.714965, - 0.648438,   0.261444],
                      [0.076308,   0.420804, - 0.903936]])
        self.assertEqual(nb0s, 3)
        np.testing.assert_allclose(gpred, gutest)

    def test_getdirs5(self):
        ndirs=10
        gutest, nb0s = difunc.get_dirs(ndirs)
        gpred = np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
        self.assertEqual(nb0s, 1)
        np.testing.assert_allclose(gpred, gutest)

if __name__ == '__main__':
    unittest.main()
