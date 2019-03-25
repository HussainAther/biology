"""
We can generalize a Turing machine to two dimensions for 2-D cellular automata
and create a Turmite. It's a portmanteau of a termite (because of how the head
of the automata moves) with Turing.

During each time step, the ant checks the color of the cell is it on. If it is black,
the ant turns to the right, changes the cell to white, and moves forward one space.
If the cell is white, the ant turns left, changes the cell to black, and moves forward.
"""

w = 75 # width
h = 52 # height
nsteps = 12000 # number of steps

# python one-line class statements
class Dir: up, right, down, left = range(4)
class Turn: left, right = False, True
class Color: white, black = ".", "#"
M = [[Color.white] * width for _ in range(height)]

x = w // 2
y = h // 2
dir = Dir.up # start upward

i = 0
while i < nsteps and 0 <= x < w and 0 <= y < h:
    turn = Turn.left if M[y][x] == Color.black else Turn.right # pythonic one-liners for turns
    M[y][x] = Color.white if M[y][x] == Color.black else Color.black
 
    dir = (4 + dir + (1 if turn else -1)) % 4
    dir = [Dir.up, Dir.right, Dir.down, Dir.left][dir]
    if   dir == Dir.up:    y -= 1
    elif dir == Dir.right:
        x -= 1
    elif dir == Dir.down:
        y += 1
    elif dir == Dir.left:
        x += 1
    else:
        assert False
    i += 1

print ("\n".join("".join(row) for row in M))


