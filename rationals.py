class RationalField():
    """The rational number field, QQ"""
    def __init__(self):
        pass

    def one(self):
        return Rational(1, 1)
    def zero(self):
        return Rational(0, 1)
    def coerce(self, x):
        # "coerces" a variable of one type into a Rational
        # TODO: possibly allow coercion from other types
        if type(x) is Rational:
            return x
        elif type(x) is int:
            return Rational(x, 1)
        else:
            raise ValueError(f"Can't coerce value {x} to Rational.")


class Rational():
    """Simple implementation of rational numbers"""
    def __init__(self, num, den):
        # input validation
        if type(num) is not int or type(den) is not int or den == 0:
            raise ValueError('num and den must both be integers with nonzero den')
        
        # internal values
        self.num = num
        self.den = den
        self.field = RationalField()

        # reduce fraction, if possible
        self._reduce()
    
    ####################
    # Internal methods #
    ####################
    
    def _reduce(self):
        g = self._gcd(self.num, self.den)

        self.num, self.den = int(self.num/g), int(self.den/g)

    def _gcd(self, a, b):
        # TODO: implement gcd algorithm instead of using the math library
        import math
        return math.gcd(a,b)
    
    def __add__(self, other):
        other = self.field.coerce(other)

        # adds this object to other and returns a new instance representing the sum
        num2, den2 = other.num, other.den

        new_num = self.num*den2 + num2*self.den
        new_den = self.den*den2

        return Rational(new_num, new_den)
    
    def __sub__(self, other):
        # subtraction is addition
        other = self.field.coerce(other)
        return self + (-1)*other
    
    def __mul__(self, other):
        # multiply self with other and return a new instance
        other = self.field.coerce(other)
        return Rational(self.num*other.num, self.den*other.den)
    
    def __rmul__(self, other):
        # sometimes I want my rational to be on the *right* but we're commutative
        return self.__mul__(other)

    def __truediv__(self, other):
        # division is multiplication
        if other == self.field.zero():
            return ZeroDivisionError()
        
        return self.__mul__(Rational(other.den, other.num))
    
    def __pow__(self, n):
        # TODO: possibly implement this for integer or (some) rational numbers
        raise NotImplementedError('Exponentiation is not implemented')
    
    def __eq__(self, other):
        # as a special case, we will allow comparison to integers
        if type(other) is int:
            other = self.field.coerce(other)

        # otherwise if they are not the right type, they are not equal
        if type(other) is not Rational:
            return False
        
        return self.num*other.den == self.den*other.num

    def __repr__(self):
        # returns a string representation of this fraction
        if self.den == 1:
            # no need to print the denominator
            return str(self.num)
        else:
            return f'{self.num}/{self.den}'
