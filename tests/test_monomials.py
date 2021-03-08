import pytest
from random import randint, seed
from groebner.monomials import MonomialOrdering, Monomial

# Can be enabled for reproducibility
# seed(12345)

# Parameters and custom testing classes
class Trash:
    def __init__(self, x):
        self.value = x
    def __repr__(self):
        return f'Trash({self.value})'

NUM_VARS = 6
CUSTOM_VARS = ['3', True, float(2), 3, Trash('can'), Trash]
CUSTOM_VAR_STR = ['3', 'True', '2.0', '3', 'Trash(can)', "<class 'test_monomials.Trash'>"]

# Tests

def grlex_lt(lst1, lst2):
    # standalone implementation of grlex order, comparing two lists of powers.

    # First compare total degree
    if sum(lst1) < sum(lst2):
        return True
    elif sum(lst1) > sum(lst2):
        return False
    else:
        # next do lexicographic order to break ties
        for i in range(len(lst1)):
            if lst1[i] < lst2[i]:
                return True
            elif lst1[i] > lst2[i]:
                return False
        # if we haven't returned yet, they are the same monomial
        return False

class TestMonomials:
    def test_input_data(self):
        # make sure we're feeding accurate data to tests
        assert len(CUSTOM_VAR_STR) == len(CUSTOM_VARS)
        assert NUM_VARS == len(CUSTOM_VARS)

    def test_ordering_instantiation(self):
        o = MonomialOrdering(num_vars = NUM_VARS)
        assert o.num_vars == NUM_VARS
        assert o.var_labels == ['x_0', 'x_1', 'x_2', 'x_3', 'x_4', 'x_5']

        o = MonomialOrdering(num_vars = NUM_VARS, labels = CUSTOM_VARS)
        assert o.num_vars == len(CUSTOM_VARS)
        assert o.var_labels == CUSTOM_VAR_STR
    
    def test_monomial_instantiation(self):
        o = MonomialOrdering(num_vars=4, labels=['x','y','z','w'])
        for _ in range(5):
            a = randint(0, 1000)
            b = randint(0, 1000)
            c = randint(0, 1000)
            d = randint(0, 1000)

            target_string = f'x^{a}y^{b}z^{c}w^{d}'
            m = Monomial([a,b,c,d], o)

            assert str(m) == target_string
            assert m.degrees == [a,b,c,d]
    
    def test_grlex_initial(self):
        """grlex ordering is the following:
        1 < z < y < x < z^2 < yz < y^2 < xz < xy < x^2 < z^3 < ...
        0   1   2   3    4    5     6    7    8     9    10"""
        o = MonomialOrdering(num_vars=3, labels=['x', 'y', 'z'], order_type='grlex')

        # check things are printing correctly
        for i, s in enumerate(['', 'z', 'y', 'x', 'z^2', 'yz', 'y^2', 'xz', 'xy', 'x^2', 'z^3']):
            assert str(o.idx_to_monomial(i)) == s
        
        # check the actual degrees are correct
        for i, d in enumerate([[0,0,0], [0,0,1], [0,1,0], [1,0,0], [0,0,2], [0,1,1], [0,2,0], [1,0,1], [1,1,0], [2,0,0], [0,0,3]]):
            assert o.idx_to_monomial(i).degrees == d
    
    def test_grlex_order(self):
        o = MonomialOrdering(num_vars=4, order_type='grlex')
        for _ in range(25):
            a = randint(0, 3000)
            b = randint(0, 3000)
            m = o.idx_to_monomial(a)
            n = o.idx_to_monomial(b)

            # check implementation of order vs reference above
            assert (m < n) == grlex_lt(m.degrees, n.degrees)
            # check that the indices have the same property
            assert (a < b) == grlex_lt(m.degrees, n.degrees)
        
        # testing equality and le
        m = o.idx_to_monomial(42)
        n = o.idx_to_monomial(43)
        p = o.idx_to_monomial(42)

        assert m != n and m <= n and n >= m and not n < m
        assert m == p and m <= p and p >= m and not p < m
    
    def test_to_idx(self):
        # assuming our idx_to_monomial method is correct, this tests to_idx()
        #   works for reasonable values
        # TODO switch to MonomialOrdering method once implemented
        o = MonomialOrdering(num_vars=4, order_type='grlex')
        for _ in range(1000):
            i = randint(0, 10000)
            m = o.idx_to_monomial(i)
            assert m.to_idx() == i
    
    def test_monomial_multiplication(self):
        """ x^4y^3z^43 * xz^2 == x^5y^3z^45 """
        o = MonomialOrdering(num_vars=3)
        m = Monomial([4, 3, 43], o)
        n = Monomial([1, 0, 2], o)
        p = Monomial([5, 3, 45], o)

        assert m*n == p
