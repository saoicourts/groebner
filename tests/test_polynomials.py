import pytest
from groebner import QQ
from groebner.polynomials import Polynomial, PolynomialRing
from groebner.rationals import Rational


class TestPolynomials:
    ORDERINGS = ['lex', 'grlex', 'grevlex']
    def test_identities(self):
        R = PolynomialRing(labels=['x','y','z'])
        r = R.random()

        assert r == r * R.one()
        assert r == R.one() * r
        assert 0 == r * R.zero()
        assert 0 == R.zero() * r
    
    def test_coercion(self):
        R = PolynomialRing(labels=['x','y','z'])
        r = R.random()

        # coerce from polynomial
        s = R.coerce(r)
        assert r == s 

        # coerce from field
        assert R.one() == R.coerce(R.field.one())

        # coerce from monomial
        m = R.ordering.random()
        p = R.coerce(m)
        assert m == p.LM()

        # coerce from numbers
        a = 42
        b = (12345)/(56789)
        p = R.coerce(a)
        with pytest.warns(RuntimeWarning):
            q = R.coerce(b)

        assert p._total_deg() == 0
        assert q._total_deg() == 0
        assert p.LC() == a
        assert (q.LC().num)/(q.LC().den) == b

        # Implicit coercions often used in practice
        assert R.one() == 1
        assert 1 == R.one()
        assert R.zero() == 0
        assert 0 == R.zero()
    
    @pytest.mark.parametrize('order', ORDERINGS)
    def test_inequality(self, order):
        R = PolynomialRing(labels=['x','y','z'], order=order)
        f = R.random()
        g = R.random()

        # Should just inherit from monomials
        assert (f <= g) == (f.LM() <= g.LM())
        # Check silly edge cases
        assert f <= f and f >= f and not f < f and not f > f
    
    def test_exponentiation(self):
        R = PolynomialRing(labels=['x','y'])
        x, y = R.get_vars()
        # semi-random example computed with sage
        f = sum([
            Rational(-8, 13)*x**4*y,
            Rational(53,54)*x**2*y**2,
            Rational(-5, 26)*y
        ])
        g = sum([
            Rational(-32768, 371293)*x**20*y**5, 
            Rational(542720, 771147)*x**18*y**6,
            Rational(-3595520, 1601613)*x**16*y**7,
            Rational(11910160, 3326427)*x**14*y**8,
            Rational(-51200, 371293)*x**16*y**5,
            Rational(-39452405, 13817466)*x**12*y**9,
            Rational(678400, 771147)*x**14*y**6,
            Rational(418195493, 459165024)*x**10*y**10,
            Rational(-1123600, 533871)*x**12*y**7,
            Rational(7443850, 3326427)*x**10*y**8,
            Rational(-32000, 371293)*x**12*y**5,
            Rational(-197262025, 221079456)*x**8*y**9,
            Rational(106000, 257049)*x**10*y**6,
            Rational(-351125, 533871)*x**8*y**7,
            Rational(18609625, 53222832)*x**6*y**8,
            Rational(-10000, 371293)*x**8*y**5,
            Rational(66250, 771147)*x**6*y**6,
            Rational(-1755625, 25625808)*x**4*y**7,
            Rational(-3125, 742586)*x**4*y**5,
            Rational(165625, 24676704)*x**2*y**6,
            Rational(-3125, 11881376)*y**5
        ])

        assert f**5 == g
    
    def test_string_representations(self):
        R = PolynomialRing(labels=['x','y'])
        x, y = R.get_vars()
        
        assert str(R.zero()) == '0'
        assert str(R.one()) == '1'
        assert str(42 * R.one()) == '42'

        assert (str(3*x**2*y-Rational(4, 3)*x*y**2)) == '3x^2y - 4/3xy^2'