import sys

"""
Create an error-free genome assembly method to create a 
circular genome from a given input list of DNA transcripts.
"""

class suffixarray:
    def __init__(self, text):
        """
        Initialize a suffix array of the text.
        """
       self.seq = self.buildsuffix(text)

    def cc(self, S, seq ):
        """
        Create the class of characters.
        """
        l = len(S)
        cclass = [0] * l
        cclass[seq[0]] = 0
        for i in range(1, l):
            if S[seq[i]] != S[seq[i - 1]]:
                cclass[seq[i]] = cclass[seq[i - 1]] + 1
            else:
                cclass[seq[i]] = cclass[seq[i - 1]]
        return cclass

    def sortit(self, S):
        """
        Determine the sorted string. 
        """
        l = len(S)
        seq = [0] * l
        count = dict()
        for i in range(l):
            count[S[i]] = count.get(S[i], 0) + 1
        charList = sorted(count.keys())
        pChar = charList[0]
        for char in charList[1:]:
            count[char] += count[pChar]
            pChar = char
        for i in range(l - 1, -1, -1):
            c = S[i]
            count[c] = count[c] - 1
            seq[count[c]] = i
        return seq

    def dsort(self, S, L, seq, classa):
        """
        Sorting algorithm to get the right order works on 
        string S, sequence seq, and class classa. 
        """
        los = len(S) # length of string
        count = [0] * los
        orderv2 = [0] * los
        for i in range(los):
            count[classa[i]] += 1
        for j in range(1, los):
            count[j] += count[j - 1]
        for i in range(los - 1, -1, -1):
            start = (seq[i] - L + los) % los
            cl = classa[start]
            count[cl] -= 1
            orderv2[count[cl]] = start
        return orderv2

    def update(self, orderv2, classa, L ):
        """
        Update the classes.
        """
        n = len(orderv2)
        n = [0] * n
        n[orderv2[0]] = 0
        for i in range(1, n):
            c = orderv2[i]
            p = orderv2[i - 1]
            mid = c + L
            midPrev = (p + L) % n
            if classa[c] != classa[p] or classa[mid] != classa[midPrev]:
                n[c] = n[p] + 1
            else:
                n[c] = n[p]
        return n

    def buildsuffix(self, S ):
        """
        Build a suffix array for the string S.
        """
        los = len(S) # length of string
        seq = self.sortit(S)
        classa = self.cc(S, seq)
        L = 1
        while L < los:
            seq = self.dsort(S, L, seq, classa)
            classa = self.update(seq, classa, L)
            L = 2 * L
        return seq

class runit:
    def __init__(self):
        reads = self.readData()
        genome = self.assembly(reads)
        print(genome)

    def readData(self):
        """
        Read the input text file of transcripts.
        """
        return list(set(sys.stdin.read().strip().split()))

    def getbwt(self, text, seq, alphabet=["$", "A", "C", "G", "T"] ):
        """
        Burrows–Wheeler transform of the suffix array.
        """
        l = len(text)
        bwt = [""] * l
        for i in range(l):
            bwt[i] = text[(seq[i] + l - 1) % l]

        counts = dict()
        starts = dict()
        for char in alphabet:
            counts[char] = [0] * (l + 1)
        for i in range(l):
            cuurent_Char = bwt[i]
            for char, count in counts.items():
                counts[char][i + 1] = counts[char][i]
            counts[cuurent_Char][i + 1] += 1
        cent_index = 0
        for char in sorted(alphabet):
            starts[char] = cent_index
            cent_index += counts[char][l]
        return bwt, starts, counts

    def getoverlap(self, text, patterns, k=12 ):
        """
        Return the longest overlap for the list of patterns
        in the text.
        """
        seq = suffixarray(text).seq
        bwt, starts, counts = self.getbwt(text, seq)
        l = len(text) - 1
        occs = dict()
        for i, pattern_list in enumerate(patterns):
            pattern = pattern_list[:k]
            top = 0
            bottom = len(bwt) - 1
            index = len(pattern) - 1
            while top <= bottom:
                if index >= 0:
                    symbol = pattern[index]
                    index -= 1
                    if counts[symbol][bottom + 1] - counts[symbol][top] > 0:
                        top = starts[symbol] + counts[symbol][top]
                        bottom = starts[symbol] + counts[symbol][bottom + 1] - 1
                    else:
                        break
                else:
                    for j in range(top, bottom + 1):
                        if not seq[j] in occs:
                            occs[seq[j]] = []
                        occs[seq[j]].append(i)
                    break
        overlap = 0
        for pos, iList in sorted(occs.items()):
            for i in iList:
                if text[pos:-1] == patterns[i][:l - pos]:
                    return i, l - pos
        return i, overlap

    def assembly(self, reads):
        """
        Assemble!
        """
        index = 0 # start here
        genome = reads[0] # output genome
        fread = reads[index] # first read
        while True:
            cead = reads[index] # cent read
            if 1 == len(reads):
                break
            del reads[index]
            index, overlap = self.getoverlap(cead + "$", reads)
            genome += reads[index][overlap:]
        index, overlap = self.getoverlap(reads[0] + "$", [fread])
        if overlap > 0:
            return genome[:-overlap]
        else:
            return genome

if __name__ == "__main__":
    runit()
