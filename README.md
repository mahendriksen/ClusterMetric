# ClusterMetric

Some notes:
So far, this 'package' consists of two programs: distancechecker and rankchecker, and one module, treedefs. All must be downloaded to the same folder.
distancechecker.py will find the distance between pairs of trees (if from a text file, all pairs consisting of one tree from file 1 and one tree from file 2)
rankchecker.py will find the rank of single trees.
Both generate a .txt file that records the number of each rank/distance in a list (so 0th entry is all rank/distance 1, 1st entry is rank/distance 2, etc.).

Both programs have several options:
1) Random (0) or Fixed (1)
	a) if 0, you have the option of generating a number of random phylogenetic trees.
	b) if 1, the program will read trees from a .txt file you provide. Details of the formatting are below
2) How many leaves in each tree?
	a) Needs to be same for all trees under comparison, and you must tell the program how many there are even if you are providing a text file
3) How many iterations? (RANDOM ONLY)
	a) How many trees do you wish to generate?
4) Binary (0) or Non-Binary (1)? (RANDOM ONLY)
	a) Generate a uniformly distributed selection of trees drawn from all trees or just binary trees. Note Non-Binary can generate binary trees, but the reverse is not true.


For the .txt files
1) You will need to write your trees in Newick form, with each tree on a new line, in one .txt file for rankchecker and two separate files for distancechecker.
2) OTU's must be encoded as integers in tuple form, e.g. (1,). For instance, (dog, (cat, fish)) and ((dog, cat), fish) must be ((1,),((2,),(3,))) and (((1,),(2,)),(3,)) respectively.
3) The only input check done is to make sure each tree has the same number of leaves. Ensure that all trees are formatted correctly as the program will not warn you.
