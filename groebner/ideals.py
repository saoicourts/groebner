from groebner.rings import RingElement
from groebner.algorithms import buchberger


class Ideal:
    """Parent class for ideals in some ring"""
    def __init__(self, ring):
        self.ring = ring
    
    def ideal_add(self, other):
        raise NotImplementedError

    def intersect(self, other):
        raise NotImplementedError
    
    def __add__(self, other):
        if not isinstance(other, Ideal):
            raise TypeError(f"Can't add types Ideal and {type(other)}.")
        if self.ring != other.ring:
            raise TypeError('Base rings for ideals must be the same to add.')
        return self.ideal_add(other)
    
class IdealFromGenerators(Ideal):
    """An ideal given by generators (elements in parent ring)"""
    def __init__(self, gens):
        # wrap in a list (if not already a list)
        if type(gens) is not list:
            gens = [gens]
        # We want elements to be instances of ring elements of the same ring
        ring = gens[0].ring
        errors = [not isinstance(x, RingElement) or x.ring != ring for x in gens]
        if sum(errors) > 0:
            raise TypeError('Input must be RingElement or list of RingElements.')

        super().__init__(ring)
        self.gens = gens

    def ideal_add(self, other):
        # if we have a list of generators, we just add concatenate them!
        if not isinstance(other, IdealFromGenerators):
            raise TypeError('Other ideal must be of type IdealFromGenerators.')
        return IdealFromGenerators(self.gens + other.gens)
    
    def intersect(self, other):
        # TODO: this
        raise NotImplementedError

class GroebnerIdeal(IdealFromGenerators):
    """An ideal given by generators that form a Groebner basis"""
    def __init__(self, gens):
        # use some algorithm to compute a Groebner basis from the given gens
        self.basis = buchberger(gens)

        super().__init__(self.basis)