import cv2
image = cv2.imread('glacier.exr', cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)

m = 256
l = 1
image2 = cv2.resize(image, (l * 2, l * 2))
for i in range(0, 3):
	for ii in range(0, l):
		for jj in range(0, l):
			image2[ii][jj][i] = 0.6 * m
			image2[ii][jj + l][i] = 0.2 * m
			image2[ii + l][jj][i] = 0.04 * m
			image2[ii + l][jj + l][i] = 0.16 * m
cv2.imwrite('ttt.exr', image2)

l2 = 1024
image3 = cv2.resize(image, (l2, l2))
cv2.imwrite('glacier-small.exr', image3)