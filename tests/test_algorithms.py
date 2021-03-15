import pytest
from groebner import QQ
from groebner.algorithms import buchberger, division_algorithm
from groebner.polynomials import PolynomialRing
from groebner.rationals import Rational



class TestAlgorithms:
    # Here we set what orderings we want to test below
    ORDERINGS = ['lex', 'grlex', 'grevlex']

    # Expected inputs and outputs for example 1

    def test_buchberger_grlex(self):
        R = PolynomialRing(labels=['x','y'], order='grlex')
        x, y = R.get_vars()

        gens = [
            x**3 - 2*x*y,
            x**2*y - 2*y**2 + x
        ]
        true_basis = [
            x**3 - 2*x*y,
            x**2*y - 2*y**2 + x,
            x**2,
            x*y,
            y**2 - Rational(1, 2)*x
        ]

        basis = buchberger(gens)

        # TODO if order is correct we can just compare lists.
        assert basis == true_basis
    
    @pytest.mark.parametrize('order', ORDERINGS)
    def test_division_by_self(self, order):
        R = PolynomialRing(labels=['x','y'], order=order)
        p = R.random(max_deg=20)

        [q], r = division_algorithm(p, [p])

        assert q == 1
        assert r == 0
    
    @pytest.mark.parametrize('order', ORDERINGS)
    def test_not_divisible(self, order):
        R = PolynomialRing(labels=['x','y'], order=order)
        x, y = R.get_vars()

        f = x**2*y - x
        g = y**2 + y

        [q], r = division_algorithm(f, [g])

        assert q == R.zero()
        assert r == f
    
    @pytest.mark.parametrize('order', ORDERINGS)
    def test_division_example_one_grlex(self, order):
        R = PolynomialRing(labels=['x','y'], order=order)
        x, y = R.get_vars()

        f = x**2*y+x*y**2+y**2
        g = x*y-1
        h = y**2-1

        [d1, d2], r = division_algorithm(f, [g, h])

        if order in ['lex', 'grevlex', 'grlex']:
            assert d1 == x+y 
            assert d2 == R.one() 
            assert r == x+y+1
        else:
            raise NotImplementedError
    
    @pytest.mark.parametrize('order', ORDERINGS)
    def test_division_example_two(self, order):
        R = PolynomialRing(labels=['x','y'], order=order)
        x, y = R.get_vars()

        f = x**7*y**2+x**3+y**2-y+1
        g = x*y*y-x
        h = x-y**3

        [d1, d2], r = division_algorithm(f, [g, h])

        if order in ['grevlex', 'grlex']:
            assert d1 == x**6
            assert d2 == 0
            assert r == x**7 + x**3 + y**2 - y + 1
        elif order == 'lex':
            # Done by hand whew!
            assert d1 == x**6 + x**5*y + x**4*y**2 + x**4 + x**3*y + x**2*y**2 + x**2 + 2*x*y + 2*y**2 + 2
            assert d2 == x**6 + x**5*y + x**4 + x**3*y + 2*x**2 + 2*x*y + 2
            assert r == 2*y**3 + y**2 - y + 1
        else:
            raise NotImplementedError
    
    @pytest.mark.parametrize('order', ORDERINGS)
    def test_random_grlex(self, order):
        R = PolynomialRing(labels=['x','y','z'], order=order)

        f = R.random(num_terms=10, max_deg=12, denominator_bound=10)
        g = R.random(num_terms=10, max_deg=10, denominator_bound=10)
        h = R.random(num_terms=10, max_deg=10, denominator_bound=10)
        k = R.random(num_terms=10, max_deg=10, denominator_bound=10)

        [d1, d2, d3], r = division_algorithm(f, [g, h, k])

        assert f == g*d1 + h*d2 + k*d3 + r