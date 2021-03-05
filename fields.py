from rings import Ring, RingElement


class Field(Ring):
    """Parent class for all fields as special cases of rings"""
    def __init__(self, name, element_type):
        super().__init__(name, element_type, is_commutative=True)

class FieldElement(RingElement):
    """Parent class for all elements of fields"""
    def __init__(self, field):
        super().__init__(field)
        self.field = field

    # on top of usual ring methods, we have a couple extras:
    def mul_inv(self):
        raise NotImplementedError

    def __truediv__(self, other):
        return self.__mul__(other.mul_inv())
    
    def __rtruediv__(self, other):
        n = self.mul_inv()
        return n.__mul__(other)

    def __pow__(self, power):
        try:
            a = int(power)
        except ValueError as err:
            raise ValueError('Powers of field elements must be integers.') from err

        if self == self.field.zero():
            if a == 0:
                raise ValueError('Undefined value 0^0.')
            return self.field.zero()
        if a == 0:
            return self.field.one()
        elif a < 0:
            # x^{-n} = (x^{-1})^n
            new = self.mul_inv()
            return new.__mul__(new.__pow__(abs(a) - 1))
        else:
            return self.__mul__(self.__pow__(a - 1))
