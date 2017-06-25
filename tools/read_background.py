import cv2
image = cv2.imread('mountains.exr', cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)

l2 = 1024
image3 = cv2.resize(image, (l2, l2))
cv2.imwrite('mountains-small.exr', image3)