from groebner.monomials import Monomial
from groebner.polynomials import Polynomial

def buchberger(gens):
    # gens is a list of polynomials
    # we want to return a list containing a Groebner basis (extending gens)
    raise NotImplementedError

def buchberger_fast(gens):
    raise NotImplementedError

def is_groebner(basis):
    raise NotImplementedError

def division_algorithm(dividend, divisors):
    """Runs the generalized division algorithm using polynomials in multiple variables"""
    # Input validation
    ring = dividend.ring

    if type(divisors) is not list:
        divs = [divisors]
    else:
        divs = divisors

    errors = [x not in ring for x in divs]
    if sum(errors) != 0:
        raise TypeError('Dividend and divisors must all come from the same ring.')

    # IVT page 65
    qs = [ring.zero()]*len(divisors)
    r = ring.zero()
    p = dividend

    while p != ring.zero():
        division_occurred = False
        for i, div in enumerate(divs):
            res = _divide_terms(p.LT(), div.LT())
            if res is not None:
                qs[i] = qs[i] + res
                p = p - res*div
                division_occurred = True
                break
        if not division_occurred:
            r = r + p.LT()
            p = p - p.LT()
    return qs, r

def _divide_terms(p, q):
    try:
        assert p.LT() == p and q.LT() == q
    except AssertionError:
        raise ValueError('Only implemented for polynomials of a single term.')

    if q == q.ring.zero():
        return None

    m, n = p.LM().degrees, q.LM().degrees
    a, b = p.LC(), q.LC()

    # pad out lists to be the same size
    while len(m) < len(n):
        m.append(p.ring.zero())
    while len(n) < len(m):
        n.append(p.ring.zero())
    
    degs = [m[i] - n[i] for i in range(len(m))]

    if all([x >= 0 for x in degs]):
        mon = Monomial(degs, p.LM().order)
        return Polynomial({mon: a/b}, p.ring)
    else:
        return None
    