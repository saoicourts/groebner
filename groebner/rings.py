class Ring:
    """Parent class/interface for all rings to inherit from"""
    def __init__(self, name, element_type, is_commutative=False):
        self.name = name
        self.element_type = element_type
        self.is_commutative = is_commutative

    # Stuff we need to implement in child classes

    def one(self):
        raise NotImplementedError

    def zero(self):
        raise NotImplementedError

    def __eq__(self, other):
        raise NotImplementedError    

    # Stuff that just works
    
    def coerce(self, x):
        # passthrough if they are the same type already
        if type(x) is type(self):
            return x
        raise NotImplementedError

    def __contains__(self, other):
        return type(other) is self.element_type
    
    def __repr__(self):
        return self.name

class RingElement:
    """Parent class for elements in some ring"""
    def __init__(self, ring):
        self.ring = ring

    # Stuff to implement

    def __add__(self, other):
        raise NotImplementedError

    def __mul__(self, other):
        raise NotImplementedError

    def __rmul__(self, other):
        # multiplication isn't necessarily commutative
        # only need to implement it if it isn't
        if self.ring.is_commutative:
            return self.__mul__(other)
        raise NotImplementedError

    def __neg__(self):
        raise NotImplementedError

    def __eq__(self, other):
        raise NotImplementedError

    # Stuff that just works

    def __pow__(self, power):
        try:
            a = int(power)
        except ValueError as err:
            raise ValueError('Powers of ring elements must be integers.') from err

        if a < 0:
            raise ValueError('Power of ring element must be nonnegative.')
        if a == 0:
            return self.ring.one()
        else:
            # compute recursively
            return self.__mul__(self.__pow__(a-1))

    def __radd__(self, other):
        # addition is always commutative
        return self.__add__(other)

    def __sub__(self, other):
        # a - b is a + (-b)
        if type(other) is not type(self):
            raise TypeError(
                f'Subtraction not defined for types {type(self)} and {type(other)}'
            )
        return self.__add__(other.__neg__())
    
    def __repr__(self):
        return f'Element of {self.ring}'
