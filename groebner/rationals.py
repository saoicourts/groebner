from random import randint
from groebner.fields import Field, FieldElement
from warnings import warn


class RationalField(Field):
    """The rational number field, QQ"""
    def __init__(self):
        super().__init__("the rational numbers", Rational)

    def one(self):
        return Rational(1, 1)

    def zero(self):
        return Rational(0, 1)
        
    def random(self, bound=10**10):
        # returns a "random" (for some definition of random) rational
        # between 0 and 1
        den = randint(2, bound)
        num = randint(1, den)
        return Rational(num, den)

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
            warn(
                "Coercing to Rational from float is inherently inexact and will "
                "be incorrect for large numerators and denominators. If you "
                "need to have certainty the values are correct, you should "
                "directly instantiate with Rational(num, den) with integers.",
                RuntimeWarning
            )
            for i in range(1, 10**7):
                if int(i*x) == i*x:
                    return Rational(int(i*x), i)
        else:
            raise ValueError(f"Can't coerce value {x} to Rational.")
    
    def __eq__(self, other):
        return type(other) is RationalField


class Rational(FieldElement):
    """Simple implementation of rational numbers"""
    def __init__(self, num, den):
        self.field = RationalField()
        super().__init__(self.field)

        # input validation
        if den == 0:
            raise ZeroDivisionError('Denominator cannot be zero.')
        try:
            assert int(num) == num and int(den) == den
            num, den = int(num), int(den)
            if int(den) < 0:
                self.num = -1*int(num)
            else:
                self.num = int(num)
            self.den = abs(int(den))
        except (ValueError, AssertionError):
            try:
                r = self.field.coerce(num/den)
                self.num = r.num
                self.den = r.den
            except ValueError:
                raise ValueError(f'Could not interpret {num} and {den}'
                                 ' as numerator and denominator of a rational.')

        # reduce fraction, if possible
        self._reduce()
    
    ####################
    # Internal methods #
    ####################

    def mul_inv(self):
        if self.den == self.field.zero():
            raise ZeroDivisionError
        return Rational(self.den, self.num)
    
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
    
    def __mul__(self, other):
        # multiply self with other and return a new instance
        try:
            other = self.field.coerce(other)
            return Rational(self.num*other.num, self.den*other.den)
        except:
            try: 
                return other.__rmul__(self)
            except:
                raise TypeError('No acceptable definition of multiplication')

    def __truediv__(self, other):
        # division is multiplication
        o = self.field.coerce(other)
        if o == self.field.zero():
            raise ZeroDivisionError
        
        return self.__mul__(o.mul_inv())
    
    def __eq__(self, other):
        try:
            other = self.field.coerce(other)
        except ValueError:
            # silently return false
            return False
        
        return self.num*other.den == self.den*other.num
    
    # QQ is totally ordered
    def __lt__(self, other):
        o = self.field.coerce(other)
        return self.num*o.den < self.den*o.num
    
    def __le__(self, other):
        return self.__lt__(other) or self.__eq__(other)
    
    def __gt__(self, other):
        o = self.field.coerce(other)
        return self.num*o.den > self.den*o.num
    
    def __ge__(self, other):
        return self.__gt__(other) or self.__eq__(other)
    
    # This works as normal
    def __abs__(self):
        return Rational(abs(self.num), self.den)
    
    # add some coercions to numeric types:
    def __int__(self):
        # the more pythonic thing would probably be to just convert blindly
        #   but I am afraid of silent lossy conversion
        self._reduce()
        if self.den == 1:
            return self.num
        else:
            raise ValueError(f'Rational {self} cannot be losslessly'
                             ' coerced to int.')
    
    def __float__(self):
        return float(self.num)/self.den

    def __repr__(self):
        # returns a string representation of this fraction
        if self.den == 1:
            # no need to print the denominator
            return str(self.num)
        else:
            return f'{self.num}/{self.den}'
