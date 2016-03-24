import numpy as np

numX = 100
numY = 52

depth = 50

dat = np.zeros((numY, numX))


depth = 50
for i in range(1, numY - 1):
    dat[i, 1:-1] = depth
    depth -= 1

for i in range(25, 31):
    for j in range(45, 56):
        dat[i][j] = 0

np.savetxt('island_depth.dat', dat, '%4d')

fid = open('island_boundry.dat', 'w')

for i in range(numX-2):
    fid.write(str(i + 2) + '\t51\t1\t1\n')

fid.close()