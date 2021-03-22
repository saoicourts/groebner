import pytest
from groebner.algorithms import buchberger , division_algorithm, is_groebner, buchberger_fast
from groebner.polynomials import PolynomialRing
from groebner.rationals import Rational


class TestAlgorithms:
    # Here we set what orderings we want to test below
    ORDERINGS = ['lex', 'grlex', 'grevlex']

    @pytest.mark.parametrize('order', ORDERINGS)
    def test_buchberger(self, order):
        try:
            R = PolynomialRing(labels=['x','y'], order=order)
            x, y = R.get_vars()

            gens = [
                x**3 - 2*x*y,
                x**2*y - 2*y**2 + x
            ]

            # the simplest (reduced) bases I have found
            if order == 'lex':
                true_basis = [x - 2*y**2, y**3]
            else:
                true_basis = [x**2, x*y, y**2 - Rational(1, 2)*x]

            basis = buchberger(gens)

            # make sure it is actually a groebner basis
            assert is_groebner(basis)

            # show they generate the same ideal
            for elm in basis:
                _, r = division_algorithm(elm, true_basis)
                assert r == 0

            for elm in true_basis:
                _, r = division_algorithm(elm, basis)
                assert r == 0
        except NotImplementedError:
            pass
    
    @pytest.mark.parametrize('order', ORDERINGS)
    def test_buchberger_variants_fixed(self, order):
        try:
            R = PolynomialRing(labels=['x','y'], order=order)
            x, y = R.get_vars()
            f = 2*x**2*y + 2*x**2 + y**2
            g = x**2*y**2 + 2*x**2 + 2*x*y**2 + 2*x*y + 1

            # get different bases with different methods
            basis1 = buchberger([f, g])
            basis2 = buchberger_fast([f, g])

            # check they are groebner bases
            assert is_groebner(basis1)
            assert is_groebner(basis2)

            # check ideal containment in a circular way
            for elm in basis1:
                _, r = division_algorithm(elm, basis2)
                assert r == 0

            for elm in basis2:
                _, r = division_algorithm(elm, basis1)
                assert r == 0
        except NotImplementedError:
            pass
    
    @pytest.mark.parametrize('order', ORDERINGS)
    def test_buchberger_variants_variable(self, order):
        try:
            R = PolynomialRing(labels=['x','y'], order=order)
            f = R.random(num_terms=5, max_deg=3, denominator_bound=3)
            g = R.random(num_terms=5, max_deg=3, denominator_bound=3)

            if f == f.ring.zero():
                f += 1
            if g == g.ring.zero():
                g += 1

            # get different bases with different methods
            basis1 = buchberger([f, g])
            basis2 = buchberger_fast([f, g])

            # check they are groebner bases
            assert is_groebner(basis1)
            assert is_groebner(basis2)

            # check ideal containment in a circular way
            for elm in basis1:
                _, r = division_algorithm(elm, basis2)
                assert r == 0

            for elm in basis2:
                _, r = division_algorithm(elm, basis1)
                assert r == 0
        except:
            NotImplementedError
    
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
