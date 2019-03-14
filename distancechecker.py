from __future__ import division
import sys
from treedefs import *
import time
from ast import literal_eval as make_tuple

def calculation_time(treeone,treetwo,ping_count,running_total,diameter): #finds distance between two trees as integer, rest of inputs are for statistical purposes.
	one_f =	dirty_calculate_f(treeone)
	two_f = dirty_calculate_f(treetwo)
	# print "First tree is ", treeone, "with a rank of ",  one_f
	# print "Second tree is ", treetwo, "with a rank of ",  two_f

	a = multihierarchy(treeone,treetwo)
	halfway_point = big_tree(a,ping_count) 
	ping_count = halfway_point[1]
	halfway_point = halfway_point[0]
	# print "Halfway point is ", halfway_point
	if len(halfway_point) > 0:
		halfway_f = clean_calculate_f(halfway_point)
		# print "Halfway point has rank of ", halfway_f
		# print "Distance between them is ", (one_f + two_f - 2*halfway_f)

	x = (one_f + two_f - 2*halfway_f)
	running_total += x
	if x > diameter:
		diameter = x
	distance_data[one_f + two_f - 2*halfway_f] += 1
	sys.stdout.write("\rCalculating iteration %i of %i . That's %f %%" % (i+1,m, (i+1)/m *100))
	sys.stdout.flush()
	return ping_count,running_total,diameter

p = input("Random (0) or Fixed (1)?")

if p == 0:
	m = input("How many iterations? ")
	n = input("How many leaves? ")
	r = input("Binary (0) or Non-Binary (1)?")
	X = range(1,n+1)
	ping_count = 0
	distance_data = [0]*(n**2-3*n+2)
	diameter = 0
	running_total = 0
	start_time = time.time()
	if r == 0:
		for i in xrange(m):
			treeone = random_binary_tree(X)
			treetwo = random_binary_tree(X)
			ping_count,running_total,diameter = calculation_time(treeone,treetwo,ping_count,running_total,diameter)
	if r == 1:
		for i in xrange(m):
			treeone = random_tree(X)
			treetwo = random_tree(X)
			ping_count,running_total,diameter = calculation_time(treeone,treetwo,ping_count,running_total,diameter)
	x = time.time() - start_time		
	print "\nTotal time is ", "--- ", x, " seconds ---" 
	print "Average time was ", x/m
	print "Average distance was ", running_total/m
	print "Diameter was ", diameter
	print "Degenerate percentage was %f %%" % (ping_count/m * 100)
if p == 1:
	ping_count = 0

	firstfile = open(raw_input("What is the first filename? "), "r") 
	first_list = firstfile.read().splitlines()
	firstfile.close()
	secondfile = open(raw_input("What is the second filename? "), "r") 
	second_list = secondfile.read().splitlines()
	secondfile.close()
	n = input("How many leaves in each tree? ")

	start_time = time.time()
	distance_data = [0]*(n**2-3*n+2)
	diameter = 0
	running_total = 0
	iterator = 0

	f = len(first_list)
	s = len(second_list)

	m = f*s

	for  i in xrange(f):
		for j in xrange(s):
			iterator += 1
			treeone = make_tuple(first_list[i])
			treetwo = make_tuple(second_list[j])
			ping_count,running_total,diameter = calculation_time(treeone,treetwo,ping_count,running_total,diameter)
			sys.stdout.write("\rCalculating iteration %i of %i . That's %f %%" % (iterator,m, iterator/m *100))
			sys.stdout.flush()

	x = time.time() - start_time		
	print "\nTotal time is ", "--- ", x, " seconds ---" 
	print "Average time was ", x/m
	print "Average distance was ", running_total/m
	print "Diameter was ", diameter
	print "Degenerate percentage was %f %%" % (ping_count/m * 100)
	
with open('%ileafdata.txt' % n, 'w+') as file:
    for number in distance_data:
        file.write("%i,\n" % number)
