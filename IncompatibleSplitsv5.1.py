#!/usr/bin/env python

import sys
import dendropy
from dendropy import treesplit
import numpy
import random
from itertools import izip_longest
import csv

'''This program will: Calculate number of conflicting bipartitions on target trees with respect to the reference tree (conflictList),
and the number of pairwise incompatible splits (pairwise_incompatibleList). It will store the bitmask for the conflicting
splits (conflictSplitList) and the newick strings for the conflicting splits on each target tree (newickConflictSplitList). If you ask, it will also generate
a null distribution of the expected mean number of conflicts and incompatible bipartitions by subsampling 10% from the full set of values
and calculating the mean 100,000 times. The distributions are saved to conflictMeans.txt and incompatibleMeans.txt. Finally, you can also generate a null distribution 
of the number of trees expected to be in conflict with each split in the reference tree by resampling.

Usage: ./IncompatibleSplitsv5.py refTreFile(newick) targetTreeList(ListOfNewickTrees) taxon_on_which_to_root_trees mean/null splitNull/none

All of the above arguments are required. If you choose to mean for the second to last argument, the program will return the mean number of conflicting bipartitions in the target
trees with respect to the reference trees. If you choose null, it will print the null distributions for the mean number of conflicts by resampling.

If you choose splitNull, the program will dump a dictionary with splits in the reference tree as keys and a list of the number of trees from each resampling that are conflict with each split to
refSplitConflictNull.txt.  IF you do not want splitNull, you need to enter none as the last argument.
'''


'''Check to make sure all the arguments have been entered on the command line'''

if len(sys.argv) != 6:
	print 'You need your reference tree (newick) as the first argument, the target tree list (newick) as the second argument, the name of the taxon on which you want to root the trees as the third, and whether you want the mean conflicts or a null distribution (mean or null) as your fourth argument and a null distribution of the number of trees expected to be in conflict with each split (splitNull or none) as your last argument. Try again. Read Usage comments.'
	sys.exit(-1)

'''Read in the trees and assign all global variables that may be modified later.'''

taxa = dendropy.TaxonSet()

refTre = dendropy.Tree.get_from_path(sys.argv[1], schema="newick", taxon_set=taxa, preserve_underscores=True)
targetTrees = dendropy.TreeList.get_from_path(sys.argv[2], schema="newick", taxon_set=taxa, preserve_underscores=True)
rootTaxon = str(sys.argv[3])

conflictList = []
conflictSplitList = []
newickConflictSplitList = []
newickConflictedReferenceSplitList = []
pairwise_incompatibleList = []
referenceSplitsConflicted = {}
refSplitDictNull = {}

resample = str(sys.argv[4])
resample = resample.lower()
splitNull = str(sys.argv[5])
splitNull = splitNull.lower()

'''Root reference tree on yli as in Hess and Goldman. Don't forget to update_splits.'''

def rootRef(reference_tree): 
  global refTre
  reference_tree.is_rooted = True
  reference_tree.update_splits()
  refRoot = reference_tree.find_node_with_taxon_label(rootTaxon)
  reference_tree.reroot_at_edge(refRoot.edge, update_splits=True)

'''Root inferred trees on yli as well. Don't forget to update_splits.'''

def rootTarget(targetTreeList):  
  global targetTrees  
  for trees in targetTreeList:
     trees.is_rooted = True
     trees.update_splits()
     rootNode = trees.find_node_with_taxon_label(rootTaxon)
     trees.reroot_at_edge(rootNode.edge, update_splits=True)

'''Calculate number of conflicting bipartitions on target trees with respect to the reference tree (conflictList),
and the number of pairwise incompatible splits (pairwise_incompatibleList). Store the bitmask for the conflicting
splits (conflictSplitList) and the newick strings for the conflicting splits on each target tree (newickConflictSplitList).'''

def conflictSplit_calc(reference_tree, targetTreeList):
  global   conflictList
  global conflictSplitList
  global newickConflictSplitList
  global newickConflictedReferenceSplitList
  global pairwise_incompatibleList
  
  for trees in targetTreeList:
    assert reference_tree.taxon_set is trees.taxon_set
    bm = reference_tree.taxon_set.all_taxa_bitmask()
    pairwise_incompatible = 0
    num_conflicting_with_reference = 0
    splits_conflicting_with_reference = []
    splits_conflicted_in_reference = []
    newick_splits_conflicting_with_reference = []
    newick_splits_conflicted_in_reference = []
    for split1 in reference_tree.split_edges:
      for split2 in trees.split_edges:
        pass
        if (not treesplit.is_compatible(split2, split1, bm)):
          pairwise_incompatible += 1
          if (split2 not in splits_conflicting_with_reference):
            splits_conflicting_with_reference.append(split2)
            num_conflicting_with_reference += 1
            newick_splits_conflicting_with_reference.append(reference_tree.taxon_set.split_as_newick_string(split2))
	  if (split1 not in splits_conflicted_in_reference):
	    splits_conflicted_in_reference.append(split1)
	    newick_splits_conflicted_in_reference.append(refTre.taxon_set.split_as_newick_string(split1))	  
    conflictList.append(num_conflicting_with_reference)  
    pairwise_incompatibleList.append(pairwise_incompatible) 
    conflictSplitList.append(splits_conflicting_with_reference)
    newickConflictSplitList.append(newick_splits_conflicting_with_reference)
    newickConflictedReferenceSplitList.append(newick_splits_conflicted_in_reference)

'''Construct a dictionary of reference tree splits as keys and the number of trees that contain splits in conflict with 
each key(split). Print splits with more than zero conflicts to stdout.'''

def refSplitsConflicted():
  global referenceSplitsConflicted
  
  for split in refTre.split_edges:
    splits = refTre.taxon_set.split_as_newick_string(split)
    referenceSplitsConflicted[splits] = 0
  
  for list in newickConflictedReferenceSplitList:
    for split in list:
      if split in referenceSplitsConflicted:
	    referenceSplitsConflicted[split] += 1
  
  for split in referenceSplitsConflicted:
    if referenceSplitsConflicted[split] > 0:
      print "Reference tree bipartition:", split
      print  "Number of target trees conflicting with this bipartition:", referenceSplitsConflicted[split]  



'''Construct a dictionary of reference tree splits as keys and the number of trees from a random sample of size N/10 that contain splits in conflict with 
each key(split). Do this 10000 times to generate a null distribution of the expected number of trees in conflict with each reference split.'''

def refSplitsConflictedNull():
  global refSplitDictNull
  bm = refTre.taxon_set.all_taxa_bitmask()
  for split in refTre.split_edges:
    if not treesplit.is_trivial_split(split, bm):
      splits = refTre.taxon_set.split_as_newick_string(split)
      refSplitDictNull[splits] = []
  for j in range(0, 10000):
    refSplitDict = {}
    randomConflicts = random.sample(newickConflictedReferenceSplitList, len(targetTrees)/10)
    for split in refTre.split_edges:
      splits = refTre.taxon_set.split_as_newick_string(split)
      refSplitDict[splits] = 0
    for list in randomConflicts:
      for split in list:
        if split in refSplitDict:
  	      refSplitDict[split] += 1
    for key in refSplitDict:
      if key in refSplitDictNull:
        refSplitDictNull[key].append(refSplitDict[key])
    
  with open('nullSplits.csv','wb') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(refSplitDictNull.keys())
    for row in izip_longest(*refSplitDictNull.values(), fillvalue=''):
      writer.writerow(row)
  file.close()

'''Generate a null distribution of the mean number of conflicting bipartitions and incompatible splits by randomly sampling N/10 from
conflictList and pairwise_incompatibleList, respectively. Write to file for plotting.'''

def nullSplits():
  meanConflicts = []
  for i in range(0,100000):
    x = random.sample(conflictList, len(targetTrees)/10)
    meanx = numpy.mean(x)
    meanConflicts.append(meanx)
  
  conflictMeans = open('conflictMeans.txt','wb')
  for meanConflict in meanConflicts:
    conflictMeans.write(str(meanConflict) + "\n")
  conflictMeans.close()

  meanIncompSplits = []
  for i in range(0,100000):
    z = random.sample(pairwise_incompatibleList, len(targetTrees)/10)
    meanz = numpy.mean(z)
    meanIncompSplits.append(meanz)

  incompatMeans = open('incompatibleMeans.txt','wb')
  for mean in meanIncompSplits:
    incompatMeans.write(str(mean) + "\n")
  incompatMeans.close()


if __name__ == "__main__":
  if resample == "mean" and splitNull == "none":
    rootRef(refTre)
    rootTarget(targetTrees)
    conflictSplit_calc(refTre,targetTrees)
    print "The mean number of conflicts:",numpy.mean(conflictList)
    print "The list of the number of conflicting bipartitions in each target tree:",conflictList
    refSplitsConflicted()
  elif resample == "null" and splitNull == "none":
    rootRef(refTre)
    rootTarget(targetTrees)
    conflictSplit_calc(refTre,targetTrees)    
    nullSplits()
  elif resample == "mean" and splitNull == "splitnull":
    rootRef(refTre)
    rootTarget(targetTrees)
    conflictSplit_calc(refTre,targetTrees)
    print "The mean number of conflicts:",numpy.mean(conflictList)
    print "The list of the number of conflicting bipartitions in each target tree:",conflictList
    refSplitsConflicted()
    refSplitsConflictedNull()
  elif resample == "null" and splitNull == "splitnull":
    rootRef(refTre)
    rootTarget(targetTrees)
    conflictSplit_calc(refTre,targetTrees)    
    nullSplits()
    refSplitsConflictedNull()
  else:
    print "you did not enter 'mean' or 'null' and 'splitNull' or 'none' as the last two arguments. look at usage instructions."












'''Below are some snippets of code used for troubleshooting and checking to make sure that everything works as expected.
All commented out with quotes at the end.

# refTre = dendropy.Tree.get_from_string('(A,((I,H,K,(L,(G,F))),(D,B,C,E,(S,R,(Q,(P,(M,(O,N))))))))',schema='newick', taxon_set=taxa)
# binaryTrees = dendropy.TreeList.get_from_path("binaryTrees.txt", schema="newick", taxon_set=taxa)
# unresolvedTrees = dendropy.TreeList.get_from_path("unresolvedTrees.txt", schema="newick", taxon_set=taxa)
# refTre.is_rooted = True
# refTre.update_splits()
# refRoot = refTre.find_node_with_taxon_label("A")
# refTre.reroot_at_edge(refRoot.edge, update_splits=True)

for trees in binaryTrees:
   trees.is_rooted = True
   trees.update_splits()
   rootNode = trees.find_node_with_taxon_label("A")
   trees.reroot_at_edge(rootNode.edge, update_splits=True)

for trees in unresolvedTrees:
   trees.is_rooted = True
   trees.update_splits()
   rootNode = trees.find_node_with_taxon_label("A")
   trees.reroot_at_edge(rootNode.edge, update_splits=True)

conflictList = []
conflictSplitList = []
newickConflictSplitList = []
pairwise_incompatibleList = []


for trees in binaryTrees:
  assert refTre.taxon_set is trees.taxon_set
  bm = refTre.taxon_set.all_taxa_bitmask()
  pairwise_incompatible = 0
  num_conflicting_with_reference = 0
  splits_conflicting_with_reference = []
  newick_splits_conflicting_with_reference = []
  for split1 in refTre.split_edges:
    for split2 in trees.split_edges:
      pass
      if (not treesplit.is_compatible(split2, split1, bm)):
        pairwise_incompatible += 1
        if (split2 not in splits_conflicting_with_reference):
          splits_conflicting_with_reference.append(split2)
          num_conflicting_with_reference += 1
          newick_splits_conflicting_with_reference.append(refTre.taxon_set.split_as_newick_string(split2))
  conflictList.append(num_conflicting_with_reference)  
  pairwise_incompatibleList.append(pairwise_incompatible) 
  conflictSplitList.append(splits_conflicting_with_reference)
  newickConflictSplitList.append(newick_splits_conflicting_with_reference)
    
       
import dendropy
from dendropy import treesplit

taxa = dendropy.TaxonSet()
#binaryTrees = dendropy.TreeList.get_from_path("binaryYeastTrees.txt", schema="newick", taxon_set=taxa)
#binaryTrees = dendropy.TreeList.get_from_path("binaryTrees.txt", schema="newick", taxon_set=taxa)
#unresolvedTrees = dendropy.TreeList.get_from_path("unresolvedYeastTrees.txt", schema="newick", taxon_set=taxa)
#unresolvedTrees = dendropy.TreeList.get_from_path("unresolvedTrees.txt", schema="newick", taxon_set=taxa)
#refTre = dendropy.Tree.get_from_path("HessRefTree.tre", schema="newick", taxon_set=taxa)

refTre = dendropy.Tree.get_from_string('(A,((I,H,K,(L,(G,F))),(D,B,C,E,(S,R,(Q,(P,(M,(O,N))))))))',schema='newick', taxon_set=taxa)

refTre.is_rooted = True
refTre.update_splits()
#refRoot = refTre.find_node_with_taxon_label("yli")
refRoot = refTre.find_node_with_taxon_label("A")
refTre.reroot_at_edge(refRoot.edge, update_splits=True)


#for trees in unresolvedTrees:
for trees in binaryTrees:
   trees.is_rooted = True
   trees.update_splits()
   #rootNode = trees.find_node_with_taxon_label("yli")
   rootNode = trees.find_node_with_taxon_label("A")
   trees.reroot_at_edge(rootNode.edge, update_splits=True)


conflictList = []

#for trees in unresolvedTrees:
for trees in binaryTrees:
   assert refTre.taxon_set is trees.taxon_set
   bm = refTre.taxon_set.all_taxa_bitmask()
   num_conflicting = 0
   for split1 in refTre.split_edges:
     for split2 in trees.split_edges:
       pass
       if (not treesplit.is_compatible(split1, split2, bm) and
         not treesplit.is_compatible(split2, split1, bm)):
         num_conflicting += 1
   conflictList.append(num_conflicting)

for trees in binaryTrees:
   assert refTre.taxon_set is trees.taxon_set
   bm = refTre.taxon_set.all_taxa_bitmask()
   num_conflicting = 0
   for split1 in refTre.split_edges:
     for split2 in trees.split_edges:
       pass
       if (not treesplit.is_compatible(split1, split2, bm) and
         not treesplit.is_compatible(split2, split1, bm)):
         num_conflicting += 1
   conflictList.append(num_conflicting)



assert refTre.taxon_set is binaryTrees[0].taxon_set
bm = refTre.taxon_set.all_taxa_bitmask()
num_conflicting = 0
for split1 in refTre.split_edges:
   for split2 in binaryTrees[0].split_edges:
     pass
     if (not treesplit.is_compatible(split1, split2, bm) and
       not treesplit.is_compatible(split2, split1, bm)):
       num_conflicting += 1
       print split1,split2
       print num_conflicting


assert refTre.taxon_set is binaryTrees[8].taxon_set
bm = refTre.taxon_set.all_taxa_bitmask()
num_conflicting = 0
for split1 in refTre.split_edges:
   for split2 in binaryTrees[8].split_edges:
     pass
     if (not treesplit.is_compatible(split2, split1, bm)):
       num_conflicting += 1
       print split1,split2
       


assert refTre.taxon_set is binaryTrees[8].taxon_set
bm = refTre.taxon_set.all_taxa_bitmask()
pairwise_incompatible = 0
num_conflicting_with_reference = 0
splits_conflicting_with_reference = []
for split1 in refTre.split_edges:
  for split2 in binaryTrees[8].split_edges:
     pass
     if (not treesplit.is_compatible(split2, split1, bm)):
       pairwise_incompatible += 1
       if (split2 not in splits_conflicting_with_reference):
         splits_conflicting_with_reference.append(split2)
         num_conflicting += 1
       #print split1,split2



for split in refTre.split_edges:
   bitmask = refTre.split_edges[split].split_bitmask
   treesplit.split_as_string(bitmask, 0)
   print(split)
   refTre.taxon_set.split_as_newick_string(split)
   
'''
