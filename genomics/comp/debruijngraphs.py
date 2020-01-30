import sys

"""
Create a genome assembly from De Bruijn graphs from input
transcripts text file.
"""

class runit:
    def __init__(self):
        self.eedges = 0 # explored edges
        self.unusednode = dict()
        self.path = []
        isBalanced = self._input()
        if not isBalanced:
            print("0")
        else:
            print("1")
            self.cycle()
            self.pprint()

    def _input(self):
        data = list(sys.stdin.read().strip().split())
        self.n, self.eedges = int(data[0]), int(data[1])
        self.unusedEdges = [[] for _ in range(self.n)]
        self.adj = [[] for _ in range(self.n)]
        self.outDeg = [0] * self.n
        self.inDeg = [0] * self.n
        self.adj_current_position = [0] * self.n
        for i in range(self.eedges):
            curFrom = int(data[2 * i + 2]) - 1
            curTo = int(data[2 * i + 3]) - 1
            self.adj[curFrom].append(curTo)
            self.outDeg[curFrom] += 1
            self.inDeg[curTo] += 1
        for i in range(self.n):
            if self.outDeg[i] != self.inDeg[i]:
                return False
        return True


    def updatePath(self, startPos):
        l = len(self.path) - 1
        self.path = self.path[startPos:l] + self.path[:startPos]
        for node, pos in self.unusednode.items():
            if pos < startPos:
                self.unusednode[node] = pos + l - startPos
            else:
                self.unusednode[node] = pos - startPos
        return

    def explore( self, s ):
        self.path.append(s)
        curPos = self.adj_current_position[s]
        curMaxPos = self.outDeg[s]
        while curPos < curMaxPos:
            self.adj_current_position[s] = curPos + 1
            if curPos + 1 < curMaxPos:
                self.unusednode[s] = len(self.path) - 1
            else:
                if s in self.unusednode:
                    del self.unusednode[s]
            v = self.adj[s][curPos]
            self.path.append(v)
            s = v
            curPos = self.adj_current_position[s]
            curMaxPos = self.outDeg[s]
            self.eedges -= 1
        return

    def cycle(self):
        self.explore(1)
        while self.eedges > 0:
            node, pos = self.unusednode.popitem()
            self.updatePath(pos)
            self.explore(node)
        return self.path

    def pprint(self):
        print(" ".join([str(node + 1) for node in self.path[:-1]]))

if __name__ == "__main__":
    runit()
