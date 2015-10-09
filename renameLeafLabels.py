#! /usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Rename leaf labels of a tree using a tab delimited text file. First column should contain exact name of existing leaf. Second column is desired name. Make sure you have unix line endings. WARNING: If you have taxon names that are solely denoted by numbers and there are support values on the tree you are annotating, you may have problems.')
parser.add_argument('-t','--table', required=True, nargs=1, help='required: specify the table of existing and desired names')
parser.add_argument('-i','--inputFile', required=True, nargs=1, help='required: specify the input tree file ex: my.bestSupport.tre')
parser.add_argument('-o','--output', required=True, nargs=1, help='required: specify the output file name ex: my.bestSupport.newLabels.tre')
args = parser.parse_args()


treeF = open(args.inputFile[0], 'r')
treeData = treeF.read()
treeF.close()

with open(args.table[0]) as taxon:
    count = 0
    for taxa in taxon:
      code = taxa.strip().split('\t')[0]
      fullname = taxa.strip().split('\t')[1]
      if code in treeData:
        treeData = treeData.replace(code,fullname)
        count += 1
    print str(count)+" codes replaced with full annotations in"+str(treeF) 

fullTree = open(args.output[0],'w')
fullTree.write(treeData)
fullTree.close
