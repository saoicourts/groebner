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
    # reference implementation of grlex order, comparing two lists of powers.

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
            a = randint(0, 100)
            b = randint(0, 100)
            c = randint(0, 100)
            d = randint(0, 100)

            target_string = ''
            for lbl, var in [('x', a), ('y', b), ('z', c), ('w', d)]:
                if var > 0:
                    if var == 1:
                        target_string += lbl
                    else:
                        target_string += f'{lbl}^{var}'

            m = Monomial([a,b,c,d], o)

            assert str(m) == target_string
            assert m.degrees == [a,b,c,d]
    
    def test_grlex_initial(self):
        """grlex ordering is the following:
        1 < z < y < x < z^2 < yz < y^2 < xz < xy < x^2 < z^3 < ...
        0   1   2   3    4    5     6    7    8     9    10"""
        o = MonomialOrdering(num_vars=3, labels=['x', 'y', 'z'], order_type='grlex')

        # check the ordering works
        degs = [[0,0,0], [0,0,1], [0,1,0], [1,0,0], [0,0,2], [0,1,1], [0,2,0], [1,0,1], [1,1,0], [2,0,0], [0,0,3]]
        for i in range(len(degs) - 1):
            m = Monomial(degs[i], o)
            n = Monomial(degs[i+1], o)
            assert m < n
    
    def test_grevlex_initial(self):
        """grevlex ordering is the following:
        1 < z < y < x < z^2 < yz < xz < y^2 < xy < x^2 < z^3 < ...
        0   1   2   3    4    5     6    7    8     9    10"""
        o = MonomialOrdering(num_vars=3, labels=['x', 'y', 'z'], order_type='grevlex')

        # check the ordering works
        degs = [[0,0,0], [0,0,1], [0,1,0], [1,0,0], [0,0,2], [0,1,1], [1,0,1], [0,2,0], [1,1,0], [2,0,0], [0,0,3]]
        for i in range(len(degs) - 1):
            m = Monomial(degs[i], o)
            n = Monomial(degs[i+1], o)
            assert m < n
    
    def test_lex_initial(self):
        """The order is pretty boring, so we're going to jump around a bit
        1 < z < z^2 < y < y^2 < x < xy^2 < x^2 < x^2z < x^2y"""
        o = MonomialOrdering(num_vars=3, labels=['x', 'y', 'z'], order_type='lex')

        degs = [[0,0,0], [0,0,1], [0,0,2], [0,1,0], [0,2,0], [1,0,0], [1,2,0], [2,0,0], [2,0,1], [2,1,0]]
        for i in range(len(degs)-1):
            m = Monomial(degs[i], o)
            n = Monomial(degs[i+1], o)
            assert m < n
    
    def test_order(self):
        # Using grlex for this example
        o = MonomialOrdering(num_vars=4, order_type='grlex')
        
        # testing equality and le
        m = Monomial([12, 32, 1], o)
        n = Monomial([12, 33, 0], o)
        p = Monomial([12, 32, 1], o)

        assert m != n and m <= n and n >= m and not n < m
        assert m == p and m <= p and p >= m and not p < m
    
    def test_monomial_multiplication(self):
        """ x^4y^3z^43 * xz^2 == x^5y^3z^45 """
        o = MonomialOrdering(num_vars=3)
        m = Monomial([4, 3, 43], o)
        n = Monomial([1, 0, 2], o)
        p = Monomial([5, 3, 45], o)

        assert m*n == p
