from groebner.rationals import RationalField as QQ
from groebner.monomials import MonomialOrdering, Monomial
from groebner.rings import Ring, RingElement


class PolynomialRing(Ring):
    """Represents a polynomial in some number of variables over a variety of
        fields. For base field use 'QQ', 'RR', 'CC', or an integer p designating
        the prime field F_p."""

    def __init__(self, num_vars=1, labels=None, base_field='QQ'):
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
            order_type='grlex'
        )

        self.vars = self.ordering.get_vars()
    
    def one(self):
        coefs = [self.field.zero()]*self.num_vars
        coefs[0] = self.field.one()
        return Polynomial(coefs, self)
    
    def zero(self):
        return Polynomial([self.field.zero()]*self.num_vars, self)

    def get_vars(self):
        # wrap monomials in polynomial wrappers
        ret = []
        for monomial in self.vars.values():
            idx = monomial.to_idx()
            coef = [self.field.zero()]*(idx + 1)
            coef[idx] = self.field.one()
            ret.append(Polynomial(coef, self))

        return ret
    
    def coerce(self, x):
        # "coerces" a variable of one type into a polynomial
        # TODO: possibly allow coercion from more types
        if type(x) is Polynomial:
            return x
        elif x in self.field:
            return Polynomial([x], self)
        elif type(x) is Monomial:
            idx = x.to_idx()
            coefs = [self.field.zero()]*idx + [self.field.one()]
            return Polynomial(coefs, self)
        elif type(x) is int or type(x) is float:
            return Polynomial([self.field.coerce(x)], self)
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
        if type(coefs) is not list:
            raise ValueError('Requires a list of coefficients in the base field')
        # s counts the number of coefficients that are not in the right base field
        s = sum([x not in parent_ring.field for x in coefs])
        if s > 0:
            raise ValueError('Coefficients not in proper field')

        self.field = parent_ring.field
        self.ring = parent_ring
        self.coefs = coefs
        self.order = parent_ring.ordering
    
    def copy(self):
        return Polynomial(self.coefs, self.ring)

    def _reduce(self):
        # Reduces the list to have minimal number of entries
        # amounts to unnecessary computation most of the time, but could save
        #   memory space eventually. At the moment I rely on it to compute
        #   the leading terms in a polynomial.
        for i, coef in reversed(list(enumerate(self.coefs))):
            if coef != self.field.zero():
                self.coefs = self.coefs[:i + 1]
                return
        # if we got here, we're looking at zero
        self.coefs = [0]
    
    def _total_deg(self):
        self._reduce()
        return len(self.coefs)
    
    def multidegree(self):
        self._reduce()
        return self.LM().degrees
    
    def LM(self):
        # Leading monomial
        self._reduce()
        return self.order.idx_to_monomial(len(self.coefs) - 1)
    
    def LC(self):
        self._reduce()
        return self.coefs[-1]
    
    def LT(self):
        self._reduce()
        coefs = [self.field.zero()]*(len(self.coefs)-1) + [self.coefs[-1]]
        return Polynomial(coefs, self.ring)

    def __eq__(self, other):
        if type(other) is not Polynomial:
            return False

        # check degrees
        if self._total_deg() != other._total_deg():
            return False
        
        # at this point they should already be reduced from the call to _deg above
        # diffs counts coefficients where they differ
        diffs = sum([self.coefs[i] != other.coefs[i] for i in range(self._total_deg())])
        
        return diffs == 0
    
    def __add__(self, other):
        summand = self.ring.coerce(other)
        
        # fill out lists to make them compatible
        length = max(len(self.coefs), len(summand.coefs))
        self.coefs += [self.field.zero()]*(length - len(self.coefs))
        summand.coefs += [self.field.zero()]*(length - len(summand.coefs))
        return Polynomial(
            [self.coefs[i] + summand.coefs[i] for i in range(length)],
            self.ring
        )
    
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
        p = self.ring.coerce(other)
        # TODO: there is probably a faster way to do this using degree
        new_lead = self.LM() * p.LM()
        coefs = [self.field.zero() for _ in range(new_lead.to_idx() + 1)]
        for i, c in enumerate(self.coefs):
            if c != self.field.zero():
                for j, d in enumerate(p.coefs):
                    if d != self.field.zero():
                        mon1 = self.order.idx_to_monomial(i)
                        mon2 = self.order.idx_to_monomial(j)
                        target_idx = (mon1 * mon2).to_idx()
                        coefs[target_idx] += c * d
        return Polynomial(coefs, self.ring)
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __repr__(self):
        s = ''
        first = True
        for i, coef in reversed(list(enumerate(self.coefs))):
            if coef != 0:
                if first:
                    if coef == -1:
                        s += '-'
                    elif coef == 1 and i == 0:
                        s += '1'
                    elif coef != 1: 
                        s += str(coef)
                    first = False
                else:
                    if coef > 0:
                        s += ' + '
                    else:
                        s += ' - '
                    
                    if abs(coef) != 1 or i == 0:
                        s += str(abs(coef))
                m = self.order.idx_to_monomial(i)
                s += str(m)
        return s
