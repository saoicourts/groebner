import pytest
from random import randint, seed
from groebner import QQ  # The field itself
from groebner.rationals import Rational  # for elements

# Can be enabled for reproducibility
# seed(12345)

# Parameters and helper methods

RATIONAL_UPPER_LIMIT = 10**6
RATIONAL_LOWER_LIMIT = -10**6

def rand_num():
    return randint(
        RATIONAL_LOWER_LIMIT,
        RATIONAL_UPPER_LIMIT
    )

def rand_den():
    return randint(
        1,
        RATIONAL_UPPER_LIMIT
    )

# Tests

class TestRationals:
    def test_float_coercion(self):
        # Floats are terrible
        # A couple of orders of magnitude more throws errors
        for _ in range(5):
            num = rand_num()
            den = rand_den()
            # Added a warning because floats are terrible
            with pytest.warns(RuntimeWarning):
                assert float(num/den) * QQ.one() == Rational(num, den)

    def test_int_coercion(self):
        for _ in range(5):
            x = rand_num()
            assert x * QQ.one() == Rational(x, 1)
    
    def test_mult_inverses(self):
        r = QQ.random()
        assert r * r.mul_inv() == QQ.one()
    
    def test_add_inverses(self):
        num = rand_num()
        den = rand_den()
        assert Rational(num, den) + Rational(-num, den) == QQ.zero()
    
    def test_multiplicative_unit(self):
        r = QQ.random()
        assert r * QQ.one() == r
    
    def test_additive_unit(self):
        r = QQ.random()
        assert r + QQ.zero() == r
    
    def test_mult_zero(self):
        r = QQ.random()
        assert r * QQ.zero() == QQ.zero()

    def test_div_zero(self):
        with pytest.raises(ZeroDivisionError):
            QQ.random() / 0
    
    def test_zero_denominator(self):
        with pytest.raises(ZeroDivisionError):
            Rational(rand_num(), 0)
