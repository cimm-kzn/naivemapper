import periodictable as pt
import operator
from itertools import count
from functools import reduce
from collections import Counter

def eratosthenes():
    """Yields the sequence of prime numbers via the Sieve of Eratosthenes."""
    D = {}  # map each composite integer to its first-found prime factor
    for q in count(2):  # q gets 2, 3, 4, 5, ... ad infinitum
        p = D.pop(q, None)
        if p is None:
            # q not a key in D, so q is prime, therefore, yield it
            yield q
            # mark q squared as not-prime (with q as first-found prime factor)
            D[q * q] = q
        else:
            # let x <- smallest (N*p)+q which wasn't yet known to be composite
            # we just learned x is composite, with p first-found prime factor,
            # since p is the first-found prime factor of q -- find and mark it
            x = p + q
            while x in D:
                x += p
            D[x] = p

class getMorgan():
    def __init__(self, isotop=False, stereo=False, element=True):
        self.__primes = tuple(x for _, x in zip(range(1000), eratosthenes()))
        self.__isotop = isotop
        self.__stereo = stereo
        self.__element = element

    def getMorgan(self, g):
        newlevels = {}
        countprime = iter(self.__primes)

        params = {n: (self.__primes[pt.elements.symbol(attr['element']).number] if self.__element else 1,
                      self.__primes[10 * attr['s_charge'] + attr['p_charge']] if self.__element else 1,
                      reduce(operator.mul,
                             (self.__primes[10 * (eattr.get('s_bond') or 0) + (eattr.get('p_bond') or 0)]
                              for eattr in g[n].values()), 1),
                      self.__primes[attr['isotop']] if self.__isotop and 'isotop' in attr else 1,
                      self.__primes[10 * (attr.get('s_stereo') or 0) + (attr.get('p_stereo') or 0)]
                      if self.__stereo else 1,
                      reduce(operator.mul,
                             (self.__primes[10 * (eattr.get('s_stereo') or 0) + (eattr.get('p_stereo') or 0)]
                              for eattr in g[n].values()), 1) if self.__stereo else 1)
                  for n, attr in g.nodes(data=True)}

        weights = {x: (newlevels.get(y) or newlevels.setdefault(y, next(countprime)))
                   for x, y in sorted(params.items(), key=operator.itemgetter(1))}

        oldnumb = numb = len(g)
        maxcount = 0
        stab = 0

        scaf = {}
        for n, m in g.edge.items():
            scaf[n] = tuple(m)

        while oldnumb >= numb and maxcount != 1 and stab < 3:
            oldnumb = numb
            neweights = {}
            countprime = iter(self.__primes)

            tmp = {}
            for n, m in scaf.items():
                """ if don't have neighbors use self weight
                """
                tmp[n] = reduce(operator.mul, (weights[x] for x in m), weights[n]**2)

            numb = len(set(tmp.values()))
            if numb == oldnumb:
                x = Counter(tmp.values())
                stab += 1
                maxcount = x[max(x)]
            else:
                stab = 0
                maxcount = 0

            weights = {x: (neweights.get(y) or neweights.setdefault(y, next(countprime)))
                       for x, y in sorted(tmp.items(), key=operator.itemgetter(1))}

        return weights