#!/usr/bin/env python3

import mms
import unittest
from mooseutils import fuzzyEqual, fuzzyAbsoluteEqual

class Test2DRCWalls(unittest.TestCase):
# class Test2DRCWalls():
    def test(self):
        labels = ['L2u', 'L2v', 'L2p']
        df1 = mms.run_spatial('2d-rc-no-slip-walls.i', 6, "--error", y_pp=labels, mpi=8)

        fig = mms.ConvergencePlot(xlabel='Element Size ($h$)', ylabel='$L_2$ Error')
        fig.plot(df1, label=labels, marker='o', markersize=8, num_fitted_points=3, slope_precision=1)
        fig.save('2d-rc-no-slip-walls.png')
        for key,value in fig.label_to_slope.items():
            print("%s, %f" % (key, value))
            self.assertTrue(fuzzyAbsoluteEqual(value, 2., .1))

if __name__ == '__main__':
    unittest.main(__name__, verbosity=2)
    # Test2DRCWalls().test()
    # Test2DRC().test()
    # Test1DRC().test()
