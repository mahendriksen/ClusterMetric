import random
import sys
import itertools
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from math import sqrt

def powerset(X):
    """Generate all subsets of X."""
    return itertools.chain.from_iterable(itertools.combinations(X, r) for r in range(len(X)+1))

def set_partitions(X):
    """
    Generate all set partitions of X.
    """
    X = list(X)
    if len(X) == 1:
        yield [X]
        return
    
    x = X[0]
    for Y in powerset(X[1:]):
        Y_set = set(Y)
        Z = [z for z in X[1:] if z not in Y_set]
        
        if len(Z) == 0:
            yield [X]
        
        else:
            for p in set_partitions(Z):
                yield [[x] + list(Y)] + p
				
def random_subset(k): #randomly makes a list of 0's and 1's so that at least one is different, intended to be used to pick a random subset
	if k == 1: 
		return [1]
	else:
		in_the_list = []
		for i in range(k):
			j = random.randint(0,1)
			in_the_list.append(j)
		x = sum(in_the_list)
		if x == 0 or x == k:
			in_the_list = random_subset(k)
		return in_the_list

def random_partition(k, iterable):
    results = [[] for i in range(k)]
    for value in iterable:
        x = random.randrange(k)
        results[x].append(value)
    fixed_results = sorted([i for i in results if len(i) > 0])
    if len(fixed_results) == 1:
        fixed_results = random_partition(k,iterable)
    return fixed_results
  
def all_trees(X):
    """
    Generate all phylogenetic trees on leaf set X.
    """
    if len(X) == 1:
        yield (X[0],)
    
    else:
        for p in set_partitions(X):
            if len(p) == 1:
                # trivial partition
                continue
            
            for subtrees in itertools.product(*[all_trees(Y) for Y in p]):
                yield tuple(sorted(subtrees))
				
def random_tree(X): #randomly picks a tree from all of RP(X), output in Newick form
	w = len(X)
	if w == 0:
		return
	elif w == 1:
		return tuple(X)
	else:
		clusterlist = []
		zz = random_partition(w,X)		
		for i in range(len(zz)):
			current_cluster = list(zz[i])
			clusterlist.append(random_tree(current_cluster))
		return tuple(sorted(clusterlist))
		
def random_binary_tree(X): #randomly picks a binary tree from all of RP(X), output in Newick form

	w = len(X)
	if w == 0:
		return
	elif w == 1:
		return tuple(X)
	else:
		clusterlist = []
		left_cluster = []
		right_cluster = []
		flipperlist = random_subset(w)
		for i in range(w):
			if flipperlist[i] == 0:
				left_cluster.append(X[i])
			if flipperlist[i] == 1:
				right_cluster.append(X[i])
		clusterlist.append(tuple(random_binary_tree(left_cluster)))
		clusterlist.append(tuple(random_binary_tree(right_cluster)))
		return tuple(clusterlist)
		
def intersection(lst1, lst2): #finds intersections between list 1 and each element of list 2, returns as list of lists
    lst3 = [list(filter(lambda x: x in lst1, sublist)) for sublist in lst2] 
    return lst3 

def dirty_calculate_f(t): #for when you want to calculate the rank and your tree is in Newick form, returns as integer
	leaves = leaf_count(t)
	tclusters = set(find_clusters(t))
	sumofclusters = 0
	for i in tclusters:
		sumofclusters += len(i)
	if (sumofclusters - len(tclusters) - (leaves-1)) < 0:
		print("ALERT")
		sys.exit()
	return (sumofclusters - len(tclusters) - (leaves-1))

def leaf_count(t): #if you have a tree in Newick form, this counts the leaves, returns as integer
	tclusters = set(find_clusters(t))
	i = 0
	for j in tclusters:
		if len(j) == 1:
			i += 1
	return i

def clean_calculate_f(t): #for when you want to calculate rank and your tree is a list of clusters, returns as integer
	sumofclusters = 0
	intest_cluster = 0
	for cluster in t:
		sumofclusters += len(cluster)
		if len(cluster) > intest_cluster:
			intest_cluster = len(cluster)
	return (sumofclusters - len(t) - (intest_cluster -1))

def maximal_clusters(clusters): #finds maximal clusters in a set of clusters , returns as list of clusters. Warning: will usually find X, need to chop out X if that's not the answer you want
	max = []
	for i in clusters:
		x = 0
		for j in clusters:
			if set(i).issubset(j) == True:
				x+=1
		if x == 1:
			max.append(i)
	return max
		
def multihierarchy(s,t): #finds the multihierarchy obtained by using the intersection algorithm on your two trees, , returns as list of clusters
	s_copy = list(set(find_clusters(s)))
	t_copy = list(set(find_clusters(t)))
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
	for i in range(len(inter_multi)):
		inter_multi[i]= tuple(inter_multi[i])
	return inter_multi

def big_tree(multi,ping_count): #Finds biggest tree that will fit into a multihierarchy, returns as tuple consisting of form (list of clusters, ping count), where ping count is used for tracking certain statistics
    treebigf = []
    biggest_count = 1	
    ping = False
    for cluster in multi: #if the multihierarchy is already a hierarchy we can save a lot of work, so this checks for that. Also removes repeated singletons and 2-clusters, since these don't change form of the hierarchy.
        a = multi.count(cluster)
        if a > 1 and len(cluster) < 3:
            multi = list(filter(lambda a: a != cluster, multi))
        if a > 1 and len(cluster) > 2:
            biggest_count = a
    if biggest_count == 1:
        treebigf = multi
    else:
        ping = True
        current_best = []
        current_best_score = 0
        checked_clusters = []
        repeating_clusters = []
        setofclusters =[]
        for cluster in multi: 
            copycount = multi.count(cluster)
            if copycount > 1 and cluster not in checked_clusters:
                checker = False
                for i in checked_clusters:
                    if cluster in i:
                        checker = True
                if not checker:        
                    checked_clusters.append(cluster)
                    entry = (cluster,copycount)
                    repeating_clusters.append(entry)
        for i in range(len(repeating_clusters)): #generates all possible combos for the hierarchy and tests which has biggest rank.
            cluster_tuple = repeating_clusters[i]
            cluster = cluster_tuple[0]
            copycount = cluster_tuple[1]
            if len(cluster) > 2:
                subclusters = list(itertools.combinations(cluster,max(1,len(cluster)-copycount)+1))
                setofclusters.append(subclusters)       
        choices = list(itertools.product(*setofclusters))
        for clustersets in choices: # building the intermediate clusters (we have made smallest one, can just add stuff back in one at a time)
            i=-1
            to_be_added = []
            multi_copy = list(set(multi))
            for cluster in clustersets:
                i += 1
                ladder = [list(cluster)]
                current_one = list(cluster)
                for x in range(len(multi_copy)): #need to do this for each specific rep'd cluster
                    nextcluster = multi_copy[x]
                    nextcluster = set(nextcluster)
                    if nextcluster.issubset(repeating_clusters[i][0]) and len(nextcluster) > 1:
                        intersect_list = [value for value in nextcluster if value in current_one]
                        multi_copy[x] = tuple(intersect_list)
                for j in repeating_clusters[i][0]:			
                    if j not in current_one:
                        current_one.append(j)
                        current_one.sort()
                        ladder.append(list(current_one))
                for adding_cluster in ladder:
                    if len(adding_cluster) > 0:
                        to_be_added.append(tuple(adding_cluster))
            for cluster in to_be_added:
                multi_copy.append(cluster)
            multi_copy = list(set(multi_copy))
            multi_copy = sorted([t for t in multi_copy if t != ()])
            test_score = clean_calculate_f(multi_copy)
            if test_score > current_best_score:
                current_best = multi_copy
                current_best_score = test_score
        treebigf = current_best
    if ping == True:
        ping_count += 1
    return treebigf,ping_count

def distancechecker(treeone,treetwo): #finds distance between two trees as integer.
    ping_count = 0
    m = leaf_count(treeone)
    n = leaf_count(treetwo)
	
    if m == n:
        X = range(1,m+1)
        one_f =	dirty_calculate_f(treeone)
        two_f = dirty_calculate_f(treetwo)

        a = multihierarchy(treeone,treetwo)
        halfway_point = big_tree(a,ping_count)
        halfway_point = halfway_point[0]
        halfway_f = clean_calculate_f(halfway_point)
    else:
        print("ERROR: Trees have different number of leaves")
        sys.exit()	
    return one_f + two_f - 2*halfway_f

def cohend(d1, d2): #finds Cohen's d for two lists of integers, returns as decimal
	# calculate the size of samples
	n1, n2 = len(d1), len(d2)
	# calculate the variance of the samples
	s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
	# calculate the pooled standard deviation
	s = sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
	# calculate the means of the samples
	u1, u2 = np.mean(d1), np.mean(d2)
	# calculate the effect size
	return (u1 - u2) / s

def leaf_set(t): #find leaf set of a tree as a list
	leaves = []
	if isinstance(t, (int, int)):
		return (t,)
	elif len(t) == 1:
		return [t[0]]
	else:
		for subtree in t:
			leaves.extend(leaf_set(subtree))
	return sorted(leaves)

def find_clusters(t): #finds the clusters of a Newick form tree
	if isinstance(t, (int, int)):
		return [(t,)]
	if len(t) == 1:
		return [t]
	else:
		this_cluster = [tuple(leaf_set(t))]
		for s in t:
			this_cluster += find_clusters(s)
		return sorted(this_cluster, key=lambda x: (-len(x),x))
		
def ncr(n, r): #finds number of combinations of r things out of n options
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer//denom
