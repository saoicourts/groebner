from groebner.ideals import IdealFromGenerators, GroebnerIdeal
from groebner.polynomials import PolynomialRing
from groebner.rationals import Rational

class TestIdeals:
    def test_containment(self):
        R = PolynomialRing(labels=['x','y'])
        x, y = R.get_vars()

        I = GroebnerIdeal([x**2+1, x*y-1])

        assert x*(x+y) in I

    def test_addition(self):
        R = PolynomialRing(labels=['x', 'y', 'z'])
        x, y, z = R.get_vars()
        f = x**2 + 2*y**2 - 3
        g = x*y - 5*x*y**5
        a = y**2*z + 10*z
        b = x*y**3

        I = GroebnerIdeal([f, g])
        J = IdealFromGenerators([a, b])

        K = J + I

        for x in [f, g, a, b]:
            assert x in K

        assert R.one() not in K
