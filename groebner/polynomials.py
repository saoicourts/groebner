from groebner.rationals import RationalField as QQ
from groebner.monomials import MonomialOrdering, Monomial
from groebner.rings import Ring, RingElement
from random import randint


class PolynomialRing(Ring):
    """Represents a polynomial in some number of variables over a variety of
        fields. For base field use 'QQ', 'RR', 'CC', or an integer p designating
        the prime field F_p."""

    def __init__(self, num_vars=1, labels=None, base_field='QQ', order='grlex'):
        # perform some input validation
        if base_field not in ['QQ','RR','CC'] and type(base_field) is not int:
            raise ValueError('Only allowable fields are QQ, RR, CC, or prime fields')
        
        # TODO: add more fields
        if base_field == 'QQ':
            self.field = QQ()
        else:
            raise NotImplementedError()
        
        # save a list of symbols we'll be using
        if labels is None:
            self.num_vars = num_vars
        else:
            errs = sum([type(x) is not str for x in labels])
            if type(labels) is list and errs == 0:
                self.num_vars = len(labels)
            else:
                raise TypeError('Parameter `vars` must be a list of strs.')

        # For now we're going to use the graded lexicographical ordering
        # https://en.wikipedia.org/wiki/Monomial_order#Graded_lexicographic_order
        self.ordering = MonomialOrdering(
            num_vars=self.num_vars,
            labels=labels,
            order_type=order
        )

        self.vars = self.ordering.get_vars()
    
    def one(self):
        return Polynomial(
            {self.ordering.constant_monomial(): self.field.one()}, 
            self
        )
    
    def zero(self):
        return Polynomial(
            {self.ordering.constant_monomial(): self.field.zero()}, 
            self
        )
    
    def random(self, num_terms=10, max_deg=20, denominator_bound=100):
        coefs = {}
        for _ in range(num_terms):
            mon = self.ordering.random(max_deg)
            coef = self.field.random(bound=denominator_bound)
            coefs[mon] = coef
        
        return Polynomial(coefs, self)

    def get_vars(self):
        # wrap monomials in polynomial wrappers
        mon_vars = self.ordering.get_vars().values()

        return [Polynomial({mon: self.field.one()}, self) for mon in mon_vars]
    
    def coerce(self, x):
        # "coerces" a variable of one type into a polynomial
        # TODO: possibly allow coercion from more types
        if type(x) is Polynomial:
            return x
        elif x in self.field:
            return Polynomial({self.ordering.constant_monomial(): x}, self)
        elif type(x) is Monomial:
            return Polynomial({x: self.field.one()}, self)
        elif type(x) is int or type(x) is float:
            return Polynomial(
                {self.ordering.constant_monomial(): self.field.coerce(x)},
                self
            )
        else:
            raise ValueError(f"Can't coerce value {x} to Polynomial.")

    def __repr__(self):
        return f'Polynomial ring over {self.field} with indeterminates {list(self.vars.values())}.'
    
    def __contains__(self, other):
        if type(other) is not Polynomial:
            return False
        return self.__eq__(other.ring)
    
    def __eq__(self, other):
        if type(other) is not PolynomialRing:
            return False
        return (self.field == other.field and self.ordering == other.ordering)


class Polynomial(RingElement):
    """Represents a polynomial in some polynomial ring"""
    def __init__(self, coefs, parent_ring):
        # input validation
        if type(parent_ring) is not PolynomialRing:
            raise ValueError('Parent ring must be an instance of PolynomialRing')
        if type(coefs) is not dict:
            raise ValueError('Requires a dict (Monomial => FieldElement')
        t = sum([y not in parent_ring.ordering for y in coefs.keys()])
        if t > 0:
            raise ValueError("Monomials don't come from the same order as parent_ring.")
        s = sum([x not in parent_ring.field for x in coefs.values()])
        if s > 0:
            raise ValueError('Coefficients not in proper field')

        # The empty dict is zero
        if coefs == {}:
            coefs = parent_ring.zero().coefs
        else:
            # remove zeros
            to_remove = []
            for mon in coefs:
                if coefs[mon] == parent_ring.field.zero():
                    to_remove.append(mon)
            for mon in to_remove:
                del coefs[mon]

        self.field = parent_ring.field
        self.ring = parent_ring
        self.coefs = coefs
        self.order = parent_ring.ordering
    
    def copy(self):
        return Polynomial(self.coefs, self.ring)
    
    def _total_deg(self):
        return self.LM().total_degree
    
    def multidegree(self):
        return self.LM().degrees
    
    def LM(self):
        # Leading monomial
        return sorted(self.coefs)[-1]
    
    def LC(self):
        return self.coefs[self.LM()]
    
    def LT(self):
        return Polynomial({self.LM(): self.LC()}, self.ring)

    def __eq__(self, other):
        try:
            o = self.ring.coerce(other)
        except ValueError as e:
            raise e

        # check monomials are the same
        try:
            assert set(self.coefs) == set(o.coefs)
        except AssertionError:
            return False
        
        # check coefficients match
        try:
            for mon in self.coefs:
                assert self.coefs[mon] == o.coefs[mon]
        except AssertionError:
            return False
        
        return True
    
    def __lt__(self, other):
        try:
            coefs = set(self.coefs.keys()).union(set(other.coefs.keys()))
            assert self.ring == other.ring

            for mon in reversed(sorted(coefs)):
                if mon not in self.coefs:
                    return True
                if mon not in other.coefs:
                    return False
                if self.coefs[mon] == other.coefs[mon]:
                    continue
                else:
                    return self.coefs[mon] < other.coefs[mon]
            return False
        except AttributeError:
            raise ValueError("Cannot compare Polynomials to non-polynomials.")
        except AssertionError:
            raise ValueError("Cannot compare polynomials in different rings.")
    
    def __gt__(self, other):
        return other.__lt__(self)
    
    def __ge__(self, other):
        return self == other or self.__gt__(other)
    
    def __le__(self, other):
        return self == other or self.__le__(other)
    
    def __add__(self, other):
        summand = self.ring.coerce(other)
        
        # get monomials
        monomials = set(self.coefs).union(set(summand.coefs))
        coefs = {}
        for mon in monomials:
            coef = self.field.zero()
            if mon in self.coefs:
                coef += self.coefs[mon]
            if mon in summand.coefs:
                coef += summand.coefs[mon]

            # Only add it in if it is nonzero.
            if coef != self.field.zero():
                coefs[mon] = coef

        return Polynomial(coefs, self.ring)
    
    def __hash__(self):
        # piggyback off strings
        return hash(self.__repr__())
    
    def __sub__(self, other):
        return self.__add__(-1*other)
    
    def __neg__(self):
        return -1*self
    
    def __pow__(self, power):
        try:
            power = int(power)
        except:
            raise TypeError(f'Power of polynomial must be integer. Got type {type(power)}')
        
        ret = self.ring.one()
        for _ in range(power):
            ret = self.__mul__(ret)
        return ret
    
    def __mul__(self, other):
        multiplicand = self.ring.coerce(other)
        # TODO: there is probably a faster way to do this
        coefs = {}
        for mon1 in self.coefs:
            for mon2 in multiplicand.coefs:
                new_mon = mon1 * mon2
                new_coef = self.coefs[mon1] * multiplicand.coefs[mon2]

                if new_mon not in coefs:
                    coefs[new_mon] = self.field.zero()

                coefs[new_mon] += new_coef
        return Polynomial(coefs, self.ring)
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __repr__(self):
        s = ''
        first = True
        for mon in reversed(sorted(self.coefs)):
            is_constant_term = mon.total_degree == 0
            coef = self.coefs[mon]
            if coef != 0:
                if first:
                    if coef == -1:
                        s += '-'
                    elif coef == 1 and is_constant_term:
                        s += '1'
                    elif coef != 1: 
                        s += str(coef)
                    first = False
                else:
                    if coef > 0:
                        s += ' + '
                    else:
                        s += ' - '
                    
                    if abs(coef) != 1 or is_constant_term:
                        s += str(abs(coef))
                s += str(mon)
        return s
