# python3

import sys

def db(c, m):
    try:
        _ = int(c)
        e = list(map(str, range(c)))
    except (ValueError, TypeError):
        e = c
        c = len(c)
    b = [0] * c * m
    s = []
    def d(i, j):
        if i > m:
            if m % j == 0:
                s.extend(b[1:j + 1])
        else:
            b[i] = b[i - j]
            d(i + 1, j)
            for j in range(b[i -j] + 1, c):
                b[i] = j
                d(i + 1, i)
    d(1, 1)
    return "".join(e[i] for i in s)
if __name__ == "__main__":
    print(db(2, int(input())))
