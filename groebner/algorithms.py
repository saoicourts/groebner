from groebner.monomials import Monomial
from groebner.polynomials import Polynomial
from itertools import permutations

def buchberger(gens):
    # gens is a list of polynomials
    # we want to return a list containing a Groebner basis (extending gens)
    G = gens.copy()
    while True:
        G_new = G.copy()
        for i, f in enumerate(G_new):
            for g in G_new[i+1:]:
                _, r = division_algorithm(s_poly(f, g), G)
                if r != f.ring.zero():
                    # make leading coefficient one
                    r = r * r.LC()**(-1)
                    if r not in G:
                        G.append(r)
        if G == G_new:
            return G

def reduced_buchberger(gens):
    # Same as above but with every new remainder we want to eliminate the parts
    #   of the previous generators that correspond to it.
    G = gens.copy()

    # First thing first: let's reduce the generators we were given
    G = reduce(G)

    while True:
        G_new = G.copy()
        for i, f in enumerate(G_new):
            for g in G_new[i+1:]:
                _, r = division_algorithm(s_poly(f, g), G)
                if r != f.ring.zero():
                    G_replacement = []
                    for g in G:
                        _, new_g = division_algorithm(g, r)
                        if new_g != f.ring.zero():
                            # make leading coefficient one
                            new_g = new_g * new_g.LC()**(-1)
                            if new_g not in G_replacement:
                                G_replacement.append(new_g)
                    G = G_replacement
                    # make leading coefficient one
                    r = r * r.LC()**(-1)
                    if r not in G:
                        G.append(r)
        if G == G_new:
            return reduce(G)

def buchberger_fast(gens):
    tuples = [(i, j) for i in range(len(gens)) for j in range(len(gens)) if j > i]
    G = gens.copy()
    while len(tuples) != 0:
        i, j = tuples[0]
        l = lcm(G[i].LM(), G[j].LM())
        p = G[i].LM()*G[j].LM()
        p = Polynomial({p: G[i].field.one()}, G[i].ring)
        if (l != p and not criterion(i, j, tuples, G)):
            _, r = division_algorithm(s_poly(G[i], G[j]), G)
            if r != r.ring.zero():
                # make leading coefficient one
                r = r * r.LC()**(-1)
                if r not in G:
                    G.append(r)
                    new_idx = [(i, len(G) - 1) for i in range(len(G) - 1)]
                    tuples = tuples + new_idx
        tuples.remove((i,j))
    return reduce(G)

def criterion(i, j, B, G):
    for k in range(len(G)):
        if k not in [i,j]:
            if (i, k) in B or (k,i) in B or (j, k) in B or (k,j) in B:
                continue
            else:
                m = G[k].LT()
                l = Polynomial(
                    {lcm(G[i].LM(), G[j].LM()): m.field.one()},
                    m.ring
                )
                d = _divide_terms(l, m)
                if d is None:
                    continue
                else:
                    # got it
                    return True
    return False


def reduce(polys):
    # Stolen graciously from Sage
    G = set(polys)
    while True:
        G_bar = set(G)
        for p in sorted(G_bar):
            G.remove(p)
            _, r = division_algorithm(p, sorted(G))
            if r != r.ring.zero():
                r = r*r.LC()**(-1)
                G.add(r)
        if G_bar == G:
            break
    return sorted(G)

def lcm(m, n):
    try:
        o = m.order
        assert n in o
        degrees = [max(m.degrees[i], n.degrees[i]) for i in range(len(m.degrees))]
        return Monomial(degrees, o)
    except AssertionError:
        raise ValueError("Parameters must be monomials must be associated with the same order.")

def s_poly(f, g):
    try:
        R = f.ring
        assert g in R
    except AssertionError:
        raise ValueError("Polynomials must come from the same ring")

    mon = lcm(f.LM(), g.LM())
    mon_poly = Polynomial({mon: R.field.one()}, R)

    return _divide_terms(mon_poly, f.LT())*f - _divide_terms(mon_poly, g.LT())*g

def is_groebner(gens):
    for i in range(len(gens) - 1):
        for j in range(i, len(gens)):
            found_division = False
            for perm in permutations(gens):
                _, r = division_algorithm(s_poly(gens[i], gens[j]), list(perm))
                if r == 0:
                    found_division = True
                    break
            if not found_division:
                return False
    return True

# Division 

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
    qs = [ring.zero()]*len(divs)
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
        raise ZeroDivisionError

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
    