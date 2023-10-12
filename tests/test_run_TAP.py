from django.test import TestCase
from tom_targets.tests.factories import SiderealTargetFactory
from mop.management.commands import run_TAP
import numpy as np

class TestTargetLoader(TestCase):

    def setUp(self):
        self.params = {'payload': '[[ 4.72338657e-03  7.83792642e-04  3.44260922e-02 -1.05776547e-03 -7.26734561e+00 -5.38086606e+00] [ 7.83792642e-04  1.99029697e-04  7.63902026e-03 -2.73329531e-04  -1.82032949e+00 -1.40172262e+00] [ 3.44260922e-02  7.63902026e-03  3.05964086e-01 -9.41870730e-03  -6.99882025e+01 -5.36284621e+01] [-1.05776547e-03 -2.73329532e-04 -9.41870730e-03  1.83896554e+03   2.15592353e+00  1.65241762e+00] [-7.26734561e+00 -1.82032949e+00 -6.99882025e+01  2.15592353e+00   1.67097527e+04  1.27342626e+04] [-5.38086606e+00 -1.40172262e+00 -5.36284621e+01  1.65241762e+00   1.27342626e+04  9.99307809e+03]]'}


    def test_load_covar_matrix(self):

        matrix = run_TAP.load_covar_matrix(self.params['payload'])

        assert(type(matrix) == type(np.zeros(1)))
        assert(matrix.shape == (6,6))