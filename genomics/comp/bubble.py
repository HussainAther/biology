# python3
import sys
import itertools

"""
De Bruijn bubble detect
"""

class b:

    def __init__(self, k, t, s):
        self.k = k
        self.t = t
        self.g = {}
        self.p = {}
        self.l = 0

        self.m = lambda v: self.g[v][1] > 1
        self.o = lambda v: len(self.g[v][0]) > 1

        self.bd_g(self.r(s))

    def r(self, s):
        a = lambda read: [read[j:j + self.k] for j in range(len(read) - self.k + 1)]
        return [u for read in s for u in a(read)]

    def bd_g(self, us):

        def ae(g, l, h):
            g.setdefault(l, [set(), 0])
            g.setdefault(h, [set(), 0])
            
            if h not in g[l][0]:
                g[l][0].add(h)
                g[h][1] += 1

        for u in us:
            l, h = u[:-1], u[1:]
            if l != h:
                ae(self.g, l, h)

    def cl(self):
        for k, v in self.g.items():
            if self.o(k):
                self.dfs(pt=[k], st=k, ce=k, dh=0)

        for _, dl in self.p.items():
            for pi in itertools.combinations(dl, r=2):
                if self.pd(pi):
                    self.l += 1

        return self.l

    def pd(self, pi):
        return len(set(pi[0]) & set(pi[1])) == 2

    def dfs(self, pt, st, ce, dh):
        if ce != st and self.m(ce):
            self.p.setdefault((st, ce), list()).append(pt[:])

        if dh == self.t:
            return

        for n in self.g[ce][0]:
            if n not in pt:
                pt.append(n)
                self.dfs(pt, st, n, dh + 1)
                pt.remove(n)

if __name__ == "__main__":
    da = sys.stdin.read().split()
    k, t, s = da[0], da[1], da[2:]
    print(b(int(k), int(t), s).cl())