def fully_connected(IMIN, IMAX):
    for i in range(IMIN, IMAX):
        for j in range(i+1, IMAX):
            print("{}\t{}".format(i, j));

def print_list(l):
    for x, y in l:
        print("{}\t{}".format(x, y))


for i in range(0, 10):
    for j in range(1, 3):
        print("{}\t{}".format(i, (i + j) % 10))

