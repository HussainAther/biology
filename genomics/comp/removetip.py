# python3
import itertools
import sys

"""
Remove tip with De Bruijn 
"""

class tr:
    def __init__(self, k, ds):
        self.k = k
        self.th = self.k
        self.g = {}
        self.paths = {}
        self.er = 0
        self.bdb(self.br(ds))

    def br(self, ds):
        ak = lambda ew: [ew[j:j + self.k] for j in range(len(ew) - self.k + 1)]
        return [m for ew in ds for m in ak(ew)]

    def bdb(self, k):

        def ag(g, e, h):
            g.setdefault(e, [set(), 0])
            g.setdefault(h, [set(), 0])
            if h not in g[e][0]:
                g[e][0].add(h)
                g[h][1] += 1
        for m in k:
            e, h = m[:-1], m[1:]
            if e != h:
                ag(self.g, e, h)

    def rt(self):
        for k, z in self.g.items():
            fa = None

            if len(z[0]) == 1 and z[1] == 0:
                fa = self.fain 
            elif len(z[0]) > 1:
                fa = self.faou 
            else : continue

            n = True
            while n:
                n = False
                for edge in z[0]:
                    if fa(edge, 0):
                        z[0].remove(edge)
                        self.er += 1
                        n = True
                        break

        return self.er

    def faou(self, c, d):
        if self.ou(c) > 1 or self.inn(c) > 1:
            return False
        
        if d == self.th:
            return False

        if self.ou(c) == 0:
            return True

        if self.faou(next(iter(self.g[c][0])), d + 1):
            self.g[c][0].pop()
            self.er += 1
            return True
        
        return False

    def fain(self, c, d):
        if d == self.th:
            return False

        if self.ou(c) == 0 or self.inn(c) > 1:
            return True
        
        if self.fain(next(iter(self.g[c][0])), d + 1):
            self.g[c][0].pop()
            self.er += 1
            return True
        
        return False

    def inn(self, z):
        return self.g[z][1]

    def ou(self, z):
        return len(self.g[z][0])

if __name__ == "__main__":
    k, ds = 15, sys.stdin.read().split()
    print(tr(k, ds).rt())