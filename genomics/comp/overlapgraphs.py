# python3
import sys

"""
Reconstruct genome assembly from transcripts 
of DNA fragments.
"""

class sarray:
    def __init__(self, letters):
        self.seq = self.bsa(letters)

    def sortit(self, x):
        l = len(x)
        seq = [0] * l
        r = {}
        for i in range(l):
            r[x[i]] = r.get(x[i], 0) + 1
        clist = sorted(r.keys())
        p = clist[0]
        for char in clist[1:]:
            r[char] += r[p]
            p = char
        for i in range(l-1, -1, -1):
            c = x[i]
            r[c] = r[c] - 1
            seq[r[c]] = i
        return seq

    def cg(self, x, seq):
        l = len(x)
        groupa = [0] * l
        groupa[seq[0]] = 0
        for i in range(1, l):
            if x[seq[i]] != x[seq[i-1]]:
                groupa[seq[i]] = groupa[seq[i-1]] + 1
            else:
                groupa[seq[i]] = groupa[seq[i-1]]
        return groupa        
    
    def dsort(self, x, e, seq, classa):
        sl = len(x)
        r = [0] * sl
        new = [0] * sl
        for i in range(sl):
            r[classa[i]] += 1
        for j in range(1, sl):
            r[j] += r[j-1]
        for i in range(sl-1, -1, -1):
            start = (seq[i]-e+sl) % sl
            cl = classa[start]
            r[cl] -= 1
            new[r[cl]] = start
        return new
    
    def up(self, new, classa, e):
        n = len(new)
        u = [0] * n
        u[new[0]] = 0
        for i in range(1, n):
            curr = new[i]
            prev = new[i-1]
            q = curr + e
            m = (prev + e) % n
            if classa[curr] != classa[prev] or classa[q] != classa[m]:
                u[curr] = u[prev] + 1
            else:
                u[curr] = u[prev]
        return u
    
    def bsa(self, x):
        sl = len(x)
        seq = self.sortit(x)
        classa = self.cg(x, seq)
        e = 1
        while e < sl:
            seq = self.dsort(x, e, seq, classa)
            classa = self.up(seq, classa, e)
            e = 2 * e
        return seq

class doit:
    def __init__(self):
        reads = self.datain()
        output = self.ass(reads)     
        print(output)   

    def datain(self):
        return list(set(sys.stdin.read().strip().split()))

    def makebwt(self, letters, seq, dnabases = ["$", "A", "C", "G", "T"]):
        l = len(letters)
        bwt = [""] * l
        for i in range(l):
            bwt[i] = letters[(seq[i]+l-1)%l]

        y = {}
        a = {}
        for char in dnabases:
            y[char] = [0] * (l + 1)
        for i in range(l):
            currChar = bwt[i]
            for char, r in y.items():
                y[char][i+1] = y[char][i]
            y[currChar][i+1] += 1
        cind = 0
        for char in sorted(dnabases):
            a[char] = cind
            cind += y[char][l]
        return bwt, a, y

    def getoverlap(self, letters, ps, k = 12):
        seq = sarray(letters).seq
        bwt, a, y = self.makebwt(letters, seq)        
        l = len(letters)-1

        f = {}
        for i, p in enumerate(ps):
            pt = p[:k]
            t = 0
            b = len(bwt) - 1
            cind = len(pt) - 1
            while t <= b:
                if cind >= 0:
                    sym = pt[cind]
                    cind -= 1
                    if y[sym][b+1] - y[sym][t] > 0:
                        t = a[sym] + y[sym][t]
                        b = a[sym] + y[sym][b+1] - 1
                    else:
                        break
                else:
                    for j in range(t, b + 1):
                        if not seq[j] in f:
                            f[seq[j]] = []
                        f[seq[j]].append(i)
                    break
        overlap = 0
        for pos, il in sorted(f.items()):
            for i in il:
                if letters[pos:-1] == ps[i][:l-pos]:
                    return i, l-pos
        return i, overlap

    def ass(self, reads):
        output = reads[0]
        currInd = 0
        firstRead = reads[currInd]
        while True:
            currRead = reads[currInd]
            if 1 == len(reads):
                break
            del reads[currInd]
            currInd, overlap = self.getoverlap(currRead+"$", reads)
            output += reads[currInd][overlap:]
        currInd, overlap = self.getoverlap(reads[0]+"$", [firstRead])
        if overlap > 0:
            return output[:-overlap]
        else:
            return output

if __name__ == "__main__":
    doit()
