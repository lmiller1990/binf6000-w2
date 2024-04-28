import guide as guide

def exercise_12():
  aln = guide.readClustalFile('exercise7.aln', guide.Protein_wGAP)
  tree = guide.runUPGMA(aln, 'poisson')
  tree.putAlignment(aln)
  tree.parsimony()
  # print("Sequence:\n\n", sequence)
  s = tree.strSequences(0, 4)
  print(s)

exercise_12()
