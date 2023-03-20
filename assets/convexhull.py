import math
import sys

EPSILON = sys.float_info.epsilon

'''
Given two points, p1 and p2,
an x coordinate, x,
and y coordinates y3 and y4,
compute and return the (x,y) coordinates
of the y intercept of the line segment p1->p2
with the line segment (x,y3)->(x,y4)
'''
def yint(p1, p2, x, y3, y4):
	x1, y1 = p1
	x2, y2 = p2
	x3 = x
	x4 = x
	px = ((x1*y2 - y1*x2) * (x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4)) / \
		 float((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4))
	py = ((x1*y2 - y1*x2)*(y3-y4) - (y1 - y2)*(x3*y4 - y3*x4)) / \
			float((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3-x4))
	return (px, py)

'''
Given three points a,b,c,
computes and returns the area defined by the triangle
a,b,c. 
Note that this area will be negative 
if a,b,c represents a clockwise sequence,
positive if it is counter-clockwise,
and zero if the points are collinear.
'''
def triangleArea(a, b, c):
	return (a[0]*b[1] - a[1]*b[0] + a[1]*c[0] \
                - a[0]*c[1] + b[0]*c[1] - c[0]*b[1]) / 2.0;

'''
Given three points a,b,c,
returns True if and only if 
a,b,c represents a clockwise sequence
(subject to floating-point precision)
'''
def cw(a, b, c):
	return triangleArea(a,b,c) < EPSILON;
'''
Given three points a,b,c,
returns True if and only if 
a,b,c represents a counter-clockwise sequence
(subject to floating-point precision)
'''
def ccw(a, b, c):
	return triangleArea(a,b,c) > EPSILON;

'''
Given three points a,b,c,
returns True if and only if 
a,b,c are collinear
(subject to floating-point precision)
'''
def collinear(a, b, c):
	return abs(triangleArea(a,b,c)) <= EPSILON

'''
Given a list of points,
sort those points in clockwise order
about their centroid.
Note: this function modifies its argument.
'''
def clockwiseSort(points):
	# get mean x coord, mean y coord
	xavg = sum(p[0] for p in points) / len(points)
	yavg = sum(p[1] for p in points) / len(points)
	angle = lambda p:  ((math.atan2(p[1] - yavg, p[0] - xavg) + 2*math.pi) % (2*math.pi))
	points.sort(key = angle)

'''
Replace the implementation of computeHull with a correct computation of the convex hull
using the divide-and-conquer algorithm
'''



'''
below are three helper functions (one to calculate a the y-level of the tangent line, and one each for finding the minimum and maximum tangents of a hull)
'''
def calc_tan(p1, p2, x):
	m = (p2[1] - p1[1]) / (p2[0]-p1[0])
	b = p1[1] - (m * p1[0])
	y_int = m*x + b
	return y_int

def find_min_tan_pair(setA, setB):

	indexA = max(range(len(setA)), key=lambda i: setA[i][0])
	indexB = min(range(len(setB)), key=lambda i: setB[i][0])

	x_mid = (setA[indexA][0] + setB[indexB][0]) // 2

	min_tan = calc_tan(setA[indexA], setB[indexB], x_mid)

	while (calc_tan(setA[indexA], setB[(indexB-1)%len(setB)], x_mid) > min_tan) or (calc_tan(setA[(indexA+1)%len(setA)], setB[indexB], x_mid) > min_tan):
		if (calc_tan(setA[indexA], setB[(indexB-1)%len(setB)], x_mid) > min_tan):
			min_tan = calc_tan(setA[indexA], setB[(indexB-1)%len(setB)], x_mid)
			indexB = (indexB-1) % len(setB)
		else:
			min_tan = calc_tan(setA[(indexA+1)%len(setA)], setB[indexB], x_mid)
			indexA = (indexA+1) % len(setA)

	return (indexA, indexB)

def find_max_tan_pair(setA, setB):

	indexA = max(range(len(setA)), key=lambda i: setA[i][0])
	indexB = min(range(len(setB)), key=lambda i: setB[i][0])

	x_mid = (setA[indexA][0] + setB[indexB][0]) // 2

	max_tan = calc_tan(setA[indexA], setB[indexB], x_mid)

	while (calc_tan(setA[indexA], setB[(indexB+1)%len(setB)], x_mid) < max_tan) or (calc_tan(setA[(indexA-1)%len(setA)], setB[indexB], x_mid) < max_tan):
		if (calc_tan(setA[indexA], setB[(indexB+1)%len(setB)], x_mid) < max_tan):
			max_tan = calc_tan(setA[indexA], setB[(indexB+1)%len(setB)], x_mid)
			indexB = (indexB+1) % len(setB)
		else:
			max_tan = calc_tan(setA[(indexA-1)%len(setA)], setB[indexB], x_mid)
			indexA = (indexA-1) % len(setA)

	return (indexA, indexB)

'''
merge takes in two sets of points, which are each already in a clockwise convex hull, finds the pair of points which are tangent to the maximum line of the hull, and minimum, then builds the new hull clockwise around those points
'''
def merge(setA, setB):
	setC = []

	max_tan_pair = find_max_tan_pair(setA, setB)
	min_tan_pair = find_min_tan_pair(setA, setB)

	indexA = min_tan_pair[0]
	indexB = max_tan_pair[1]

	setC.append(setA[indexA])
	while indexA != max_tan_pair[0]:
		indexA = (indexA + 1) % len(setA)
		setC.append(setA[indexA])
	setC.append(setB[indexB])
	while indexB != min_tan_pair[1]:
		indexB = (indexB + 1) % len(setB)
		setC.append(setB[indexB])

	return setC


'''
This is the base case function - used on smallest number of points after dividing full set

'''

def baseCase(points):
	hull = []
	for p1 in points:
		for p2 in points:
			if p1 == p2:
				break
			valid = True
			for p3 in points:
				if p1 == p3 or p2 == p3:
					break
				elif collinear(p1, p2, p3):
					continue
				elif ccw(p1, p2, p3):
					valid = False
					break
			if valid:
				if p1 not in hull:
					hull.append(p1)
				if p2 not in hull:
					hull.append(p2)

	clockwiseSort(hull)
	return hull


'''
integer k represents the constant at which to use the base case for merge
'''
def mergeHull(points, k):
	if len(points) <= k:
		return baseCase(points)
	
	midpoint = len(points) // 2
	left = mergeHull(points[:midpoint], k)
	right = mergeHull(points[midpoint:], k)

	return merge(left, right)


def computeHull(points):
	x = lambda p: (p[0])
	points.sort(key = x)

	return mergeHull(points, 5)
