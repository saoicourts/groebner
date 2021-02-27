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
        # TODO: possibly allow coercion from more types
        if type(x) is Rational:
            return x
        elif type(x) is int:
            return Rational(x, 1)
        elif type(x) is float:
            # This is a bit of a lazy way to do this
            # Only works with denominators less than 10^7
            # TODO: Make this better?
            for i in range(1, 10**7):
                if int(i*x) == i*x:
                    return Rational(int(i*x), i)
        else:
            raise ValueError(f"Can't coerce value {x} to Rational.")
    
    def __contains__(self, item):
        if type(item) is not Rational:
            return False
        else:
            return self.__eq__(item.field)
    
    def __eq__(self, other):
        return type(other) is RationalField
    
    def __repr__(self):
        return "the rational numbers"


class Rational():
    """Simple implementation of rational numbers"""
    def __init__(self, num, den):
        self.field = RationalField()

        # input validation
        if den == 0:
            raise ValueError('Denominator cannot be zero.')
        if int(num) != num or int(den) != den:
            try:
                r = self.field.coerce(num/den)
                self.num = r.num
                self.den = r.den
            except ValueError:
                raise ValueError(f'Could not interpret {num} and {den}'
                                    ' as numerator and denominator of a rational.')
        else:
            if int(den) < 0:
                self.num = -1*int(num)
            else:
                self.num = int(num)
            self.den = abs(int(den))

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
    
    def __radd__(self, other):
        return self.__add__(other)
    
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
        o = self.field.coerce(other)
        if o == self.field.zero():
            raise ZeroDivisionError()
        
        return self.__mul__(Rational(o.den, o.num))
    
    def __rtruediv__(self, other):
        recip = Rational(self.den, self.num)
        return recip * other
    
    def __pow__(self, n):
        # TODO: possibly implement this for (some) rational powers
        if type(n) is Rational and n.den == 1:
            n  = n.num

        if type(n) is int:
            power = abs(n)
            if n < 0:
                out = Rational(self.den, self.num)
            else:
                out = Rational(self.num, self.den)
        else:
            raise NotImplementedError('Exponentiation is not implemented for'
                                      f' type {type(n)}.')

        for _ in range(power - 1):
            out = out.__mul__(out)
        return out
    
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
