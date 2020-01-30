# python3

import sys

class eup:
    def __init__(self, adj):
        self.adj = adj
        self.ua()

    def ua(self):
        self.n = len(self.adj)
        self.na = 0 
        self.no = {} 
        self.id = {}
        self.od = {}
        self.ac = {}
        self.path = []
        self.un = []
        for w, vl in self.adj.items():
            self.id[w] = self.id.get(w, 0)
            for v in vl:
                self.id[v] = self.id.get(v, 0) + 1
            l = len(vl)
            self.od[w] = l
            self.na += l
            self.ac[w] = 0
    
    def ae(self):
        if type(self.adj) is dict:
            for v in self.adj.keys():
                if self.id[v] != self.od[v]:
                    if self.id[v] < self.od[v]:
                        self.un.append(v)
                    else:
                        self.un.insert(0, v)
            if len(self.un) > 0:
                self.adj[self.un[0]].append(self.un[1])
                self.od[self.un[0]] += 1
                self.id[self.un[1]] += 1
            return    
        for v in range(self.n):
            if self.id[v] != self.od[v]:
                if self.id[v] < self.od[v]:
                    self.un.append(v)
                else:
                    self.un.insert(0, v)
        if len(self.un) > 0:
            self.adj[self.un[0]].append(self.un[1])
            self.od[self.un[0]] += 1
            self.id[self.un[1]] += 1
        return
    
    def explore(self, s):
        self.path.append(s)
        cu = self.ac[s]
        cmp = self.od[s]
        while cu < cmp:
            self.ac[s] = cu + 1
            if cu + 1 < cmp:
                self.no[s] = len(self.path) - 1
            else:
                if s in self.no:
                    del self.no[s]
            v = self.adj[s][cu]
            self.path.append(v)
            s = v
            cu = self.ac[s]
            cmp = self.od[s]
            self.na -= 1
        return

    def up(self, st):
        l = len(self.path) - 1
        self.path = self.path[st:l] + self.path[:st]
        for nd, p in self.no.items():
            if p < st:
                self.no[nd] = p + l - st
            else:
                self.no[nd] = p - st
        return

    def ec(self):
        if type(self.adj) is dict:
            w, vl = self.adj.popitem()
            self.adj[w] = vl
            self.explore(w)
        else:
            self.explore(0)
        while self.na > 0:
            nd, p = self.no.popitem()
            self.up(p)
            self.explore(nd)
        return self.path
    
    def calculateeup(self):
        self.ae()
        self.ec()
        if len(self.un) > 0:
            for i in range(len(self.path)-1):
                if self.path[i] == self.un[0] and self.path[i+1] == self.un[1]:
                    self.up(i+1)
                    break
        return self.path          

class StringReconstruction:
    def __init__(self):
        self.k, self.adj = self.rd()
        self.path = eup(self.adj).ec()
        print(self.rfp(self.path)[:-self.k+1])

    def rd(self):
        data = list(sys.stdin.read().strip().split())
        adj = self.db(len(data[0]), data) 
        return len(data[0]), adj
    
    def db(self, k, patterns):
        aj = {}
        for p in patterns:
            if p[:k-1] in aj:
                aj[p[:k-1]].append(p[1:])
            else:
                aj[p[:k-1]] = []
                aj[p[:k-1]].append(p[1:])
            if p[1:] not in aj:
                aj[p[1:]] = []
        return aj

    def rfp(self, path):
        return path[0] + "".join(z[-1] for z in path[1:])

if __name__ == "__main__":
    StringReconstruction()