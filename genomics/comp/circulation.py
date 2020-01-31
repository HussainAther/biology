# python3
import queue
import sys

"""
Find circulation in a network.
"""

class Edge:
    def __init__(self, u, v, l, c):
        self.u = u
        self.v = v
        self.l = l
        self.c = c
        self.diff = c - l
        self.f = 0

class fgr:

    def __init__(self, n):
        self.e = []
        self.gr = [[] for _ in range(n+2)]
        self.dv = [0] * (n+2)
        self.D = 0

    def aed(self, from_, to, l, c):
        fe = Edge(from_, to, l, c)
        be = Edge(to, from_, 0, 0)
        self.gr[from_].append(len(self.e))
        self.e.append(fe)
        self.gr[to].append(len(self.e))
        self.e.append(be)
        self.dv[from_] += l
        self.dv[to] -= l

    def size(self):
        return len(self.gr)

    def getid(self, from_):
        return self.gr[from_]

    def gedge(self, id):
        return self.e[id]

    def af(self, id, f):
        self.e[id].f += f
        self.e[id ^ 1].f -= f
        self.e[id].diff -= f
        self.e[id ^ 1].diff += f

class cp:
    def __init__(self):
        gr, n, m = self.inp()
        f, fs = self.fc(gr, n, m)
        self.pp(f, fs)
    
    def inp(self):
        vertex_count, edge_count = map(int, input().split())
        gr = fgr(vertex_count)
        for _ in range(edge_count):
            u, v, l, c = map(int, input().split())
            gr.aed(u-1, v-1, l, c)
        
        for v in range(vertex_count):
            if gr.dv[v] < 0:
                gr.aed(vertex_count, v, 0, -gr.dv[v])
            if gr.dv[v] > 0:
                gr.aed(v, vertex_count+1, 0, gr.dv[v])
                gr.D += gr.dv[v]

        return gr, vertex_count, edge_count

    def bfs(self, gr, from_, to):
        X = float("inf")
        hp = False
        n = gr.size()
        dist = [float("inf")] * n
        path = []
        parent = [(None, None)] * n
        q = queue.Queue()
        dist[from_] = 0
        q.put(from_)
        while not q.empty():
            cfn = q.get()
            for id in gr.getid(cfn):
                cre = gr.gedge(id)
                if float("inf") == dist[cre.v] and cre.diff > 0:
                    dist[cre.v] = dist[cfn] + 1
                    parent[cre.v] = (cfn, id)
                    q.put(cre.v)
                    if cre.v == to:
                        while True:
                            path.insert(0, id)
                            cx = gr.gedge(id).diff
                            X = min(cx, X)
                            if cfn == from_:
                                break
                            cfn, id = parent[cfn]
                        hp = True
                        return hp, path, X
        return hp, path, X

    def mf(self, gr, from_, to):
        f = 0
        while True:
            hp, path, X = self.bfs(gr, from_, to)
            if not hp:
                return f
            for id in path:
                gr.af(id, X)
            f += X
        return f
    
    def fc(self, gr, n, m):
        f = self.mf(gr, n, n+1)
        fs = [0] * m
        if f != gr.D:
            return False, fs
        else:
            for i in range(m):
                currForwardEdge = gr.e[i*2]
                fs[i] = currForwardEdge.f + currForwardEdge.l
            return True, fs

    def pp(self, f, fs):
        if not f:
            print("NO")
        else:
            print("YES")
            print("\n".join(map(str, fs)))

if __name__ == "__main__":
    cp()
