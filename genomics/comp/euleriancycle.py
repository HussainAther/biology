# python3
import sys

"""
Find an Eulerian cycle (circuit) from input graph data of a directed graph.
"""

class EulerianCycle:
    def __init__(self):
        self.nue = 0 
        self.nwue = dict() 
        self.path = []
        isb = self.inp()
        if not isb:
            print("0")
        else:
            print("1")    
            self.ec()
            self.pr()

    def inp(self):
        data = list(sys.stdin.read().strip().split())
        self.n, self.nue = int(data[0]), int(data[1])
        self.adj = [[] for _ in range(self.n)]
        self.ue = [[] for _ in range(self.n)]
        self.od = [0] * self.n
        self.id = [0] * self.n
        self.acp = [0] * self.n
        for i in range(self.nue):
            cf = int(data[2*i+2])-1
            ct = int(data[2*i+3])-1
            self.adj[cf].append(ct)
            self.od[cf] += 1
            self.id[ct] += 1
        for i in range(self.n):
            if self.od[i] != self.id[i]:
                return False
        return True
    
    def e(self, s):
        self.path.append(s)
        cp = self.acp[s]
        cmp = self.od[s]
        while cp < cmp:
            self.acp[s] = cp + 1
            if cp + 1 < cmp:
                self.nwue[s] = len(self.path) - 1
            else:
                if s in self.nwue:
                    del self.nwue[s]
            v = self.adj[s][cp]
            self.path.append(v)
            s = v
            cp = self.acp[s]
            cmp = self.od[s]
            self.nue -= 1
        return

    def up(self, sp):
        l = len(self.path) - 1
        self.path = self.path[sp:l] + self.path[:sp]
        for n, p in self.nwue.items():
            if p < sp:
                self.nwue[n] = p + l - sp
            else:
                self.nwue[n] = p - sp
        return

    def ec(self):
        self.e(1)
        while self.nue > 0:
            n, p = self.nwue.popitem()
            self.up(p)
            self.e(n)
        return self.path

    def pr(self):
        print(" ".join([str(n+1) for n in self.path[:-1]]))       

if __name__ == "__main__":
    EulerianCycle()
