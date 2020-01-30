# python3
import sys

class sa:
    def __init__(self, text):
        self.seq = self.bsa(text)

    def sortio(self, v):
        l = len(v)
        seq = [0] * l
        co = {}
        for i in range(l):
            co[v[i]] = co.get(v[i], 0) + 1
        clist = sorted(co.keys())
        p = clist[0]
        for char in clist[1:]:
            co[char] += co[p]
            p = char
        for i in range(l-1, -1, -1):
            c = v[i]
            co[c] = co[c] - 1
            seq[co[c]] = i
        return seq

    def cc(self, v, seq):
        l = len(v)
        c = [0] * l
        c[seq[0]] = 0
        for i in range(1, l):
            if v[seq[i]] != v[seq[i-1]]:
                c[seq[i]] = c[seq[i-1]] + 1
            else:
                c[seq[i]] = c[seq[i-1]]
        return c        
    
    def dsort(self, v, e, seq, groupa):
        sl = len(v)
        co = [0] * sl
        new = [0] * sl
        for i in range(sl):
            co[groupa[i]] += 1
        for j in range(1, sl):
            co[j] += co[j-1]
        for i in range(sl-1, -1, -1):
            start = (seq[i]-e+sl) % sl
            cl = groupa[start]
            co[cl] -= 1
            new[co[cl]] = start
        return new
    
    def up(self, new, groupa, e):
        n = len(new)
        newClass = [0] * n
        newClass[new[0]] = 0
        for i in range(1, n):
            curr = new[i]
            prev = new[i-1]
            mid = curr + e
            midPrev = (prev + e) % n
            if groupa[curr] != groupa[prev] or groupa[mid] != groupa[midPrev]:
                newClass[curr] = newClass[prev] + 1
            else:
                newClass[curr] = newClass[prev]
        return newClass
    
    def bsa(self, v):
        sl = len(v)
        seq = self.sortio(v)
        groupa = self.cc(v, seq)
        e = 1
        while e < sl:
            seq = self.dsort(v, e, seq, groupa)
            groupa = self.up(seq, groupa, e)
            e = 2 * e
        return seq

def im(a1, a2, eps = 5):
    n = 0
    for i in range(len(a1)):
        if a1[i] != a2[i]:
            n += 1
            if n > eps:
                return False
    return True

class doitagain:
    def __init__(self):
        r = self.rd()
        output = self.ass(r)     
        print(output)   

    def rd(self):
        data = list(sys.stdin.read().strip().split())
        td = set()
        for i in range(1, len(data)):
            if im(data[i], data[i-1]):
                td.add(i)
        r = [d for i, d in enumerate(data) if not i in td]
        return r

    def makebwt(self, text, seq, dnabases = ["$", "A", "C", "G", "T"]):
        l = len(text)
        bwt = [""] * l
        for i in range(l):
            bwt[i] = text[(seq[i]+l-1)%l]

        cos = {}
        b = {}
        for char in dnabases:
            cos[char] = [0] * (l + 1)
        for i in range(l):
            cr = bwt[i]
            for char, co in cos.items():
                cos[char][i+1] = cos[char][i]
            cos[cr][i+1] += 1
        cindex = 0
        for char in sorted(dnabases):
            b[char] = cindex
            cindex += cos[char][l]
        return bwt, b, cos

    def flo(self, text, pat, k = 5):
        seq = sa(text).seq
        bwt, b, cos = self.makebwt(text, seq)        
        l = len(text)-1

        o = {}
        for i, p in enumerate(pat):
            pt = p[:k]
            t = 0
            b = len(bwt) - 1
            cindex = len(pt) - 1
            while t <= b:
                if cindex >= 0:
                    sy = pt[cindex]
                    cindex -= 1
                    if cos[sy][b+1] - cos[sy][t] > 0:
                        t = b[sy] + cos[sy][t]
                        b = b[sy] + cos[sy][b+1] - 1
                    else:
                        break
                else:
                    for j in range(t, b + 1):
                        if not seq[j] in o:
                            o[seq[j]] = []
                        o[seq[j]].append(i)
                    break
        ov = 0
        for pos, il in sorted(o.items()):
            for i in il:
                if im(text[pos:-1], pat[i][:l-pos]):
                    return i, l-pos
        return i, ov

    def ass(self, r):
        output = r[0]
        cind = 0
        fread = r[cind]
        while True:
            cu = r[cind]
            if 1 == len(r):
                break
            del r[cind]
            cind, ov = self.flo(cu+"$", r)
            output += r[cind][ov:]
        cind, ov = self.flo(r[0]+"$", [fread])
        if ov > 0:
            return output[:-ov]
        else:
            return output

if __name__ == "__main__":
    doitagain()