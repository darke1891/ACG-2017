import pyopenvdb as vdb
import numpy
filename = '75.vdb'

names = ['density', 'flame', 'fuel', 'heat', 'react']

for name in names:
	grid = vdb.read(filename, name)
	data = {}
	for iter in grid.iterOnValues():
		if not iter.active:
			print (iter.active, iter.depth, iter.count, iter.value)
		data[iter.min] = iter.value
	(range0, range1) = grid.evalActiveVoxelBoundingBox()
	print (name, len(data.keys()))
	print (max(data.values()), min(data.values()))
	print (range0, range1)
	fout = open(name + '.voxel', 'w')
	for i in range(0, 3):
		fout.write(str(range1[i]) + ' ')
	fout.write('\n')
	for i in range(1, range1[0] + 1):
		for j in range(1, range1[1] + 1):
			for k in range(1, range1[2] + 1):
				index = (i, j, k)
				value = data.pop(index, 0)
				fout.write(str(value) + ' ')
			fout.write('\n')
	fout.close()
