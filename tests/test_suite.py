import pytest
from groebner import QQ  # The field itself
from groebner.rationals import Rational  # for elements

class TestRationals:
    def test_float_coercion(self):
        # Check that we coerce floats correctly
        assert (142/18373)*QQ.one() == Rational(142, 18373)
    def test_int_coercion(self):
        # Check the same for ints
        assert 40 * Rational(2, 17) == Rational(80, 17)
    