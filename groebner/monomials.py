from math import factorial
from warnings import warn


class MonomialOrdering():
    """Has information pertaining to the ordering of monomials as well as 
        conversions between lists of coefficients and polynomials"""
    def __init__(self, num_vars, labels=None, order_type='grlex'):
        self.num_vars = num_vars
        
        if labels is None:
            # if no variable labels are provided, we make our own
            self.var_labels = []
            for i in range(num_vars):
                self.var_labels.append(f'x_{i}')
        else:
            # if the user does something silly here it might not look great,
            #   but shouldn't break things.
            if type(labels) is not list:
                self.var_labels = [str(labels)]
            else:
                self.var_labels = list(map(str, labels))
        
        if order_type == 'grlex':
            self.idx_to_monomial = self._idx_to_mon_grlex
            self.monomial_to_idx = self._mon_to_idx_grlex
        else:
            raise NotImplementedError('Only implemented for "grlex" ordering.')
        
        # TODO Implement grevlex, which should be faster.

        # always assume variables have decreasing order
        # x_1 > x_2 > x_3 > ... > x_n
        # In total degree k, there are (n+k-1) C (n-1) monomials

    # we want to expose the "variables" (monomials of degree 1)
    def get_vars(self):
        vars = {}
        for i in range(self.num_vars):
            lst = [0]*(self.num_vars)
            lst[i] = 1
            vars[self.var_labels[i]] = Monomial(lst, self)
        return vars
    
    def _mon_to_idx_grlex(self, mon):
        #TODO 
        pass

    def _idx_to_mon_grlex(self, idx):
        try:
            idx = int(idx)
        except:
            raise TypeError('Index must be castable to integer value')

        if idx == 0:
            return Monomial([0]*self.num_vars, self)
        
        degree, remainder = self._get_total_degree(idx)
        # remainder is between 1 and (n+degree-1)C(n-1), inclusive

        # start off with the "lowest" state and at each step do the thing that 
        # raises the energy the least
        base = [0]*(self.num_vars - 1) + [degree]
        for _ in range(int(remainder) - 1):
            # find lowest nonzero power
            lowest_idx = self.num_vars - 1
            while base[lowest_idx] == 0:
                lowest_idx -= 1

            # take one from this variable and add one to the next largest
            base[lowest_idx] -= 1
            base[lowest_idx - 1] += 1

            # the rest of the elements from the lowest index go to the last var
            base[-1] = base[lowest_idx]
            if lowest_idx != len(base) - 1:
                base[lowest_idx] = 0
        
        return Monomial(base, self)
        
    def _get_total_degree(self, idx):
        # "spin off" graded pieces (total degree)
        if idx == 0:
            return 0, 0
        i = 1
        while idx > self._choose(self.num_vars + i - 1, self.num_vars - 1):
            idx -= self._choose(self.num_vars + i - 1, self.num_vars - 1)
            i += 1
        
        return i, idx
    
    def __eq__(self, other):
        if type(other) is not MonomialOrdering:
            return False
        
        return (self.num_vars == other.num_vars and
                self.var_labels == other.var_labels)
    
    def _choose(self, n, k):
        return factorial(n)/(factorial(k)*factorial(n-k))

class Monomial():
    """Wrapper for a list that represents a monomial"""
    def __init__(self, degrees, order):
        is_not_int = lambda x: type(x) is not int
        errors = sum(list(map(is_not_int, degrees)))
        if type(degrees) is not list or errors > 0:
            raise TypeError("Monomial only accepts a list of ints.")

        self.degrees = degrees
        self.total_degree = sum(degrees)
        self.order = order
        self.num_vars = self.order.num_vars
        # TODO self.to_idx = self.order.monomial_to_idx
    
    def to_idx(self):
        # TODO remove this once implemented in the ordering
        # automatically validates the input
        base = self.degrees[:]
        total_degree = sum(base)

        # This will be the same process but reverse
        remainder = 0
        while base[-1] < total_degree:
            remainder += 1
            # find lowest nonzero power (besides the last variable)
            lowest_idx = self.num_vars - 2
            while base[lowest_idx] == 0:
                lowest_idx -= 1

            # take one from this variable and add one to the next largest
            base[lowest_idx] -= 1
            if self.num_vars - 1 != lowest_idx + 1:
                base[lowest_idx + 1] += 1 + base[-1]
                base[-1] = 0
            else:
                base[-1] += 1

        idx = 0
        for i in range(total_degree):
            idx += self._choose(self.num_vars + i - 1, self.num_vars - 1)
        
        return int(idx + remainder)

    def __mul__(self, other):
        # for now we only allow multiplication with other monomials
        if type(other) is not Monomial:
            raise TypeError('Monomials can only be multipled by other monomials.')
        if self.order != other.order:
            raise ValueError('Can only multiply monomials with compatible'
                             'monomial orderings.')

        degrees = [self.degrees[i] + other.degrees[i] for i in range(self.num_vars)]

        return Monomial(degrees, self.order)
    
    def __eq__(self, other):
        if type(other) is not Monomial:
            return False
        return self.degrees == other.degrees and self.order == other.order
    
    def __lt__(self, other):
        if type(other) is not Monomial:
            raise TypeError('Can only compare monomials with others.')
        return self.to_idx() < other.to_idx()
    
    def __le__(self, other):
        return self.__eq__(other) or self.__lt__(other)
    
    def __gt__(self, other):
        if type(other) is not Monomial:
            raise TypeError('Can only compare monomials with others.')
        return other.__lt__(self)
    
    def __ge__(self, other):
        return self.__eq__(other) or self.__gt__(other)
    
    def __repr__(self):
        s = ''
        var = self.order.get_vars()
        for i, lbl in enumerate(var.keys()):
            if self.degrees[i] > 0:
                s += lbl
                if self.degrees[i] > 1:
                    s += '^' + str(self.degrees[i])
        return s

    def _choose(self, n, k):
        return factorial(n)/(factorial(k)*factorial(n-k))