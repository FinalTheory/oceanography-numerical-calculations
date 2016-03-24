numX = 100
numY = 52

depth = 50

fid = open('rectangle_depth.dat', 'w')

for i in range(numY):
    if i == 0 or i == numY - 1:
        fid.write('0   ' * numX + '\n')
    else:
        fid.write('0' + ('%4d' % depth) * ( numX - 2) + '   0' + '\n')
        depth -= 1

fid.close()

fid = open('rectangle_boundry.dat', 'w')

for i in range(numX-2):
    fid.write(str(i + 2) + '\t51\t1\t1\n')

fid.close()