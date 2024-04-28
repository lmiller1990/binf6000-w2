import guide as guide

# # s1 = Sequence('THISLINE-', Protein_wGAP)
# # s2 = Sequence('ISALIGNED', Protein_wGAP)
# # aln = Alignment([s1, s2])

b62 = guide.readSubstMatrix('Files_WS2/blosum62.matrix', guide.Protein_Alphabet)

# # print(aln)
# seqs = readFastaFile('dao.fa', Protein_Alphabet)
# last = seqs[-1]

# print(b62)
# aln1 = align(seqs[-1], seqs[0], b62, -7)
# aln1.writeHTML("play.html")

s1 = "AAAAAAA"
s2 = "AAAAACC"
s3 = "ACCCCCA"
s4 = "CCCCC"
s5 = "AAACCCC"

# guide.align(s1, s2, )

seqs = map(lambda x: guide.Sequence(x, guide.DNA_Alphabet, gappy=True),[s1,s2,s3,s4])
seq1 = guide.Sequence(s1, guide.DNA_Alphabet, gappy=True, name="s1")
seq2 = guide.Sequence(s2, guide.DNA_Alphabet, gappy=True, name="s2")
seq3 = guide.Sequence(s3, guide.DNA_Alphabet, gappy=True, name="s3")
seq4 = guide.Sequence(s4, guide.DNA_Alphabet, gappy=True, name="s4")
seq5 = guide.Sequence(s5, guide.DNA_Alphabet, gappy=True, name="s5")


def calc_score(aln: guide.Alignment) -> int:
  s = 0
  s1 = aln.seqs[0]
  s2 = aln.seqs[1]
  for i in range(len(aln.seqs[0])):
    if s1[i] == "-" or s2[i] == "-":
      s -= 1
    elif s1[i] == s2[i]:
      s += 1
  return s

def handle(s1: guide.Sequence, s2: guide.Sequence):
  aln = guide.align(s1, s2, b62)
  print(f"Score: {calc_score(aln)}.\n{str(aln).strip()}\n")

pairs = [
  [seq1, seq2],
  [seq1, seq3],
  [seq1, seq4],
  [seq1, seq5],

  [seq2, seq3],
  [seq2, seq4],
  [seq2, seq5],

  [seq3, seq4],
  [seq3, seq5],

  [seq4, seq5],
]

for pair in pairs:
  s1, s2 = pair
  handle(s1, s2)


# c = guide.readClustalFile("./wip.out", guide.DNA_Alphabet)
# o = guide.runUPGMA(c, "poisson")
# print(o)

