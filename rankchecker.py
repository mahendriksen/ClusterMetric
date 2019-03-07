from __future__ import division
import sys
from treedefs import *
from ast import literal_eval as make_tuple
import time

p = input("Random (0) or Fixed (1)?")

n = input("How many leaves in each tree? ")
X = range(1,n+1)

distance_data = [0]*int(((n-1)*(n-2)/2) +1)

running_total = 0


if p==0:
	m = input("How many iterations? ")
	r = input("Binary (0) or Non-Binary (1)?")
	start_time = time.time()
	X = range(1,n+1)
	if r == 0:
		for s in xrange(m):
			treeone = random_binary_tree(X)
			one_f =	dirty_calculate_f(treeone)
			running_total += one_f
			distance_data[one_f] += 1
			sys.stdout.write("\rCalculating iteration %i of %i . That's %f %%" % (s+1,m, (s+1)/m *100))
			sys.stdout.flush()
	if r == 1:
		for s in xrange(m):
			treeone = random_tree(X)
			one_f =	dirty_calculate_f(treeone)
			running_total += one_f
			distance_data[one_f] += 1
			sys.stdout.write("\rCalculating iteration %i of %i . That's %f %%" % (s+1,m, (s+1)/m *100))
			sys.stdout.flush()
	with open('%ileafranks.txt' % n, 'w+') as file:
		for number in distance_data:
			file.write("%i,\n" % number)

if p == 1:
	filename = raw_input("What is the filename? ")
	firstfile = open(filename, "r") 
	first_list = firstfile.read().splitlines()
	firstfile.close()
	start_time = time.time()
	m = len(first_list)
	for i in xrange(m):
		treeone = make_tuple(first_list[i])
		print treeone
		one_f =	dirty_calculate_f(treeone)
		print one_f
		running_total += one_f
		distance_data[one_f] += 1
		sys.stdout.write("\rCalculating iteration %i of %i . That's %f %%" % (i+1,m, (i+1)/m *100))
		sys.stdout.flush()
	filename = filename[:-4]
	with open('%sdata.txt' % filename, 'w+') as file:
		for number in distance_data:
			file.write("%i,\n" % number)
		
x = time.time() - start_time		
print "Total time is ", "--- ", x, " seconds ---" 
print "Average time was ", x/m
print "Average rank was ", running_total/m

print distance_data

# t = ((((4,), ((10,), (1,))), (9,)), ((7,), (((8,), ((3,), (2,))), ((5,), (6,)))))

# tclusters = set(find_clusters(t))
# print tclusters
# sumofclusters = 0
# for i in tclusters:
	# sumofclusters += len(i)
# print sumofclusters - len(tclusters)