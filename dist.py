from __future__ import division
from treedefs import *
from ast import literal_eval as make_tuple
import time
import numpy

def multicheck(s,t): #finds the multihierarchy obtained by using the intersection algorithm on your two trees, , returns as list of clusters
	s_copy = s[:]
	t_copy = t[:]
	inter_multi = []
	max1 = maximal_clusters(s_copy)
	max2 = maximal_clusters(t_copy)
	while len(max1) > 0:
		for i in max1:
			inter = intersection(i,max2)
			for j in inter:
				if len(j) >0:
					inter_multi.append(j)
		for i in max1:
			s_copy.remove(i)
		for j in max2:
			t_copy.remove(j)
		max1 = maximal_clusters(s_copy)
		max2 = maximal_clusters(t_copy)
	for i in xrange(len(inter_multi)):
		inter_multi[i]= tuple(inter_multi[i])
	return inter_multi

def binding(t,cluster1,cluster2): #binds two clusters together in a tree. Checking whether a binding is valid must be done prior to this.
	t_copy = t[:]
	new_cluster = tuple(sorted(cluster1+cluster2))
	if len(cluster1) > 1:
		t_copy.remove(cluster1)
	if len(cluster2) > 1:
		t_copy.remove(cluster2)
	t_copy.append(new_cluster)
	return sorted(t_copy, key=lambda x: (-len(x),x))

def mod_max(clusters): #finds maximal clusters in a set of clusters , returns as list of clusters. Will never return just one cluster unless tree is one leaf.
	max = []
	if len(clusters) == 1:
		max = clusters
	else:
		cluster_copy = clusters[:]
		for i in clusters:
			x = 0
			for j in clusters:
				if set(i).issubset(j) == True:
					x+=1
			if x == 1:
				max.append(i)
		if len(max) == 1:
			cluster_copy.remove(max[0])
			max = mod_max(cluster_copy)
	return max
	
def inter(list1,list2): #intersection of two lists
	intersect = [a for a in list1 if a in list2]
	return intersect

def restriction(clusters,Y): #restricting clusters to a subset of their leaves
	restr = list(set([tuple(inter(Y,i)) for i in clusters if inter(Y,i) is not None if len(inter(Y,i)) != 0]))
	return restr
	
def binding_finder(clusters): #finds potential clusters for binding
	binding_tuples = []
	if len(clusters) == 1:
		return
	else:
		maxes = mod_max(clusters)
		if len(maxes) > 2:
			for i in xrange(len(maxes)):
				for j in xrange(i+1,len(maxes)):
					binding_tuples += [(maxes[i],maxes[j])]
		for i in xrange(len(maxes)):
			Y = list(maxes[i])
			subtree = restriction(clusters,Y)
			binding_tuples += binding_finder(subtree) if binding_finder(subtree) is not None else []
		return binding_tuples
				
def binding_maker(clusters): #finds all trees that are potential bindings
	binding_list = []
	binding_tuples = binding_finder(clusters)
	for i in binding_tuples:
		binding_list += [binding(clusters,i[0],i[1])]
	return binding_list

def unbinding(t,cluster1,list1,list2): #unbinds a cluster into two clusters. Must check if unbinding is valid before doing this
	t_copy = t[:]
	new_cluster1 = tuple(sorted([i for sub in list1 for i in sub]))
	new_cluster2 = tuple(sorted([i for sub in list2 for i in sub]))
	t_copy.remove(tuple(cluster1))
	t_copy.append(new_cluster1)
	t_copy.append(new_cluster2)
	t_copy = list(set(t_copy))
	return sorted(t_copy, key=lambda x: (-len(x),x))
		
def unbinding_finder(clusters): #finds all potential unbindings
	unbinding_tuples = []
	if len(clusters) == 1:
		return
	else:
		maxes = mod_max(clusters)
		for z in xrange(len(maxes)):
			Y = list(maxes[z])
			subtree = restriction(clusters,Y)
			our_maxes = mod_max(subtree)
			binary_number = '1'*len(our_maxes)
			num_its = (int(binary_number,2)+1) // 2
			for i in xrange(1,num_its):
				list1 = []
				list2 = []
				this_partition = [int(x) for x in list('0'*(len(our_maxes)-len(bin(i)[2:]))+bin(i)[2:])]
				for j in xrange(len(this_partition)):
					if this_partition[j] == 0:
						list1.append(our_maxes[j])
					if this_partition[j] == 1:
						list2.append(our_maxes[j])
				if len(list1) > 1 or len(list1[0]) == 1:
					if len(list2) > 1 or len(list2[0]) == 1:
						unbinding_tuples += [(maxes[z],list1,list2)]
			unbinding_tuples += unbinding_finder(subtree) if unbinding_finder(subtree) is not None else []
		return unbinding_tuples

def unbinding_maker(clusters): #finds all trees that can be made from unbindings
	unbinding_list = []
	unbinding_tuples = unbinding_finder(clusters)
	for i in unbinding_tuples:
		unbinding_list += [unbinding(clusters,i[0],i[1],i[2])]
	return unbinding_list
	
def calculation_time(treeone,treetwo,ping_count,running_bound_total,bound_diameter): #finds upper bound of distance between two trees as integer, rest of inputs are for statistical purposes.
	one_f =	dirty_calculate_f(treeone)
	two_f = dirty_calculate_f(treetwo)
	print "First tree is ", treeone, "with a rank of ",  one_f
	print "Second tree is ", treetwo, "with a rank of ",  two_f
	a = multihierarchy(treeone,treetwo)
	halfway_point = big_tree(a,ping_count) 
	ping_count = halfway_point[1]
	halfway_point = halfway_point[0]
		# print "Halfway point is ", halfway_point
	if len(halfway_point) > 0:
		halfway_f = clean_calculate_f(halfway_point)
		# print "Halfway point has rank of ", halfway_f
		x = (one_f + two_f - 2*halfway_f)
		print "Upper bound is ", x


	running_bound_total += x
	if x > bound_diameter:
		bound_diameter = x
	upper_bound_data[one_f + two_f - 2*halfway_f] += 1
	# sys.stdout.write("\rCalculating iteration %i of %i . That's %f %%" % (i+1,m, (i+1)/m *100))
	# sys.stdout.flush()
	return ping_count,running_bound_total,bound_diameter, x

def path_gen(clusters1,clusters2,k): #finds list of potential sequences of ups and downs (in binary) to get from one tree to the other
	path_list = []
	rank_one = clean_calculate_f(clusters1)
	rank_two = clean_calculate_f(clusters2)
	for j in xrange(k,abs(rank_one - rank_two)+2,-2):
		for i in xrange(2**j):
			number_of_ones = sum([int(q) for q in list(bin(i)[2:])])
			if 2*number_of_ones - j == rank_two - rank_one:
				if (bin(i)[2:2+j-number_of_ones]) != '0'*(j-number_of_ones):
					path_list.append('0'*(j-len(bin(i)[2:]))+bin(i)[2:])
	return path_list

def neigh_via_path(clusters1,paths): #finds trees in a given path
	tree_paths ={tuple(clusters1):''}
	paths_so_far = []
	for a in xrange(len(paths)):
		print "Part %d of %d" % (a+1,len(paths))
		for b in xrange(0,len(paths[a])):
			new_neighbours = []
			prefix = paths[a][:b]
			if paths[a][:b+1] not in paths_so_far:
				if paths[a][b] == '0':
					for x in tree_paths:
						if tree_paths[x] == prefix:
							new_neighbours += unbinding_maker(list(x))
				else:
					for x in tree_paths:
						if tree_paths[x] == prefix:
							new_neighbours += binding_maker(list(x))	
			for new_friend in new_neighbours:
				new_friend = tuple(new_friend)
				if not tree_paths.has_key(new_friend):
					tree_paths[tuple(new_friend)] = paths[a][:(b+1)]
	return tree_paths

def distance_via_paths(clusters1,clusters2,upper): #finds shortest distance among all possible paths
	shortest_distance = upper
	poss_paths = path_gen(clusters1,clusters2,upper)
	first_path = []
	second_path = []
	for given_path in poss_paths:
		a = len(given_path)
		first_path.append(''.join(given_path[:a//2]))
		second_path.append(bin(2**len(given_path[a//2:]) -int(''.join(given_path[a//2:]),base=2))[2:])
	first_path = list(set(first_path))
	second_path = list(set(first_path))
	print "Tree 1 computing"
	trees_in_first_path = neigh_via_path(clusters1,first_path)
	print "Tree 2 computing"
	trees_in_second_path = neigh_via_path(clusters2,second_path)
	for tree_neighbour in trees_in_first_path:
				if trees_in_second_path.has_key(tree_neighbour):
					if len(trees_in_first_path[tree_neighbour]) + len(trees_in_second_path[tree_neighbour]) < shortest_distance:
						shortest_distance = len(trees_in_first_path[tree_neighbour]) + len(trees_in_second_path[tree_neighbour])
	return shortest_distance
	
	
p = input("Random (0) or Fixed (1)?")

if p == 0:
	m = input("How many iterations? ")
	n = input("How many leaves? ")
	r = input("Binary (0) or Non-Binary (1)?")
	X = range(1,n+1)
	ping_count = 0
	upper_bound_data = [0]*(n**2-3*n+2)
	distance_data = [0]*(n**2-3*n+2)
	equality_counter = [0]*(n**2-3*n+2)
	percentage_data = [0]*(n**2-3*n+2)
	bound_diameter = 0
	distance_diameter = 0
	running_distance_total = 0
	running_bound_total = 0
	biggest_difference = 0
	times = []
	
	if r == 0:
		for i in xrange(m):
			print "Testing Pair ", i+1
			start_time = time.time()
			treeone = random_binary_tree(X)
			treetwo = random_binary_tree(X)
			clusters1 = find_clusters(treeone)
			clusters2 = find_clusters(treetwo)
			ping_count,running_bound_total,bound_diameter, upper = calculation_time(treeone,treetwo,ping_count,running_bound_total,bound_diameter)
			this_distance = distance_via_paths(clusters1,clusters2,upper)
			distance_data[this_distance] += 1
			running_distance_total += this_distance
			if this_distance > distance_diameter:
				distance_diameter = this_distance
			if upper - this_distance > biggest_difference:
				biggest_difference = upper - this_distance
			if upper == this_distance:
				equality_counter[upper] += 1
			print "The upper bound was ", upper, "but the true distance is ", this_distance
			x = time.time() - start_time	
			times.append(x)
	if r == 1:
		for i in xrange(m):
			print "Testing Pair ", i+1
			start_time = time.time()
			treeone = random_tree(X)
			treetwo = random_tree(X)
			clusters1 = find_clusters(treeone)
			clusters2 = find_clusters(treetwo)
			ping_count,running_bound_total,bound_diameter, upper = calculation_time(treeone,treetwo,ping_count,running_bound_total,bound_diameter)
			this_distance = distance_via_paths(clusters1,clusters2,upper)
			distance_data[this_distance] += 1
			running_distance_total += this_distance
			if this_distance > distance_diameter:
				distance_diameter = this_distance
			if upper - this_distance > biggest_difference:
				biggest_difference = upper - this_distance
			if upper == this_distance:
				equality_counter[upper] += 1
			print "The upper bound was ", upper, "but the true distance is ", this_distance
			x = time.time() - start_time	
			times.append(x)
	print "\nTotal time is ", "--- ", sum(times), " seconds ---" 
	print "Mean time was ", numpy.mean(times), " seconds ---"
	print "Median time was", numpy.median(times), " seconds ---"
	print "Average true distance was ", running_distance_total/m
	print "True diameter was ", distance_diameter
	print "Average bound distance was ", running_bound_total/m
	print "Bound diameter was ", bound_diameter
	print "Biggest difference was ", biggest_difference
	
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
			clusters1 = find_clusters(treeone)
			clusters2 = find_clusters(treetwo)
			ping_count,running_total,diameter, upper = calculation_time(treeone,treetwo,ping_count,running_total,diameter)
			this_distance = distance_via_paths(clusters1,clusters2,upper)
			print "\nThe upper bound was ", upper, "but the true distance is ", this_distance
			# sys.stdout.write("\rCalculating iteration %i of %i . That's %f %%" % (iterator,m, iterator/m *100))
			# sys.stdout.flush()

	x = time.time() - start_time		
	print "\nTotal time is ", "--- ", x, " seconds ---" 
	# print "Average time was ", x/m
	# print "Average distance was ", running_total/m
	# print "Diameter was ", diameter
	# print "Degenerate percentage was %f %%" % (ping_count/m * 100)
	
for i in xrange(len(upper_bound_data)):
	if upper_bound_data[i] > 0:
		percentage_data[i] = str(equality_counter[i]/upper_bound_data[i])
		print int(equality_counter[i])/int(upper_bound_data[i])
		


with open('%ileafboundsmore.txt' % n, 'w+') as file:
    for number in upper_bound_data:
        file.write("%i,\n" % number)
		
with open('%ileafdistancesmore.txt' % n, 'w+') as file:
    for number in distance_data:
        file.write("%i,\n" % number)
		
with open('%ileafpercentages.txt' % n, 'w+') as file:
    for number in percentage_data:
        file.write("%s\n" % number)
