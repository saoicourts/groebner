from rationals import RationalField as QQ

class PolynomialRing():
    """Represents a polynomial in some number of variables over a variety of
        fields. For base field use 'QQ', 'RR', 'CC', or an integer p designating
        the prime field F_p."""

    def __init__(self, num_vars = 1, base_field='QQ'):
        # perform some input validation
        if base_field not in ['QQ','RR','CC'] or type(base_field) is not int:
            raise ValueError('Only allowable fields are QQ, RR, CC, or prime fields')
        
        # TODO: add more fields
        if base_field == 'QQ':
            self.field = QQ()
        else:
            raise NotImplementedError()
        
        # save a list of symbols we'll be using
        self.vars = vars
    
    def one(self):
        return self.field.one()
    
    def zero(self):
        return self.field.zero()

    # internal methods

    def __repr__(self):
        return f''


class Polynomial():
    """Represents a polynomial in some polynomial ring"""
    def __init__(self, coefs, parent_ring):
        # input validation
        if type(parent_ring) is not PolynomialRing:
            raise ValueError('Parent ring must be an instance of PolynomialRing')
        if type(coefs) is not list:
            raise ValueError('Requires a list of coefficients in the base field')
        # s counts the number of coefficients that are not in the right base field
        s = sum([type(x) != parent_ring.field for x in coefs])
        if s > 0:
            raise ValueError('Coefficients not in proper field')

        self.field = parent_ring.field
        self.ring = parent_ring
        self.coefs = coefs

    # internal methods

    def _reduce(self):
        """Reduces the list to have deg(f) entries"""
        for i, coef in reversed(list(enumerate(self.coefs))):
            if coef != self.field.zero():
                self.coefs = self.coefs[:i + 1]
                return
        # if we got here, we're looking at zero
        self.coefs = [0]
    
    def _deg(self):
        self._reduce()
        return len(self.coefs)

    def __eq__(self, other):
        # first check degrees
        if self._deg() != other._deg():
            return False
        
        # at this point they should already be reduced from the call to _deg above
        # diffs counts coefficients where they differ
        diffs = sum([self.coefs[i] != other.coefs[i] for i in range(self._deg())])
        
        return diffs == 0
    
    def __add__(self, other):
        return Polynomial(
            [self.coefs[i] + other.coefs[i] for i in range(self._deg())],
            self.ring
        )
