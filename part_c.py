import guide as guide
import os

# tree = guide.readNewickFile('Files_WS2/cyp1a1.nwk')

# aln = guide.readClustalFile('Files_WS2/cyp1a1.aln', guide.Protein_Alphabet)
# tree.putAlignment(aln)
# tree.parsimony()

# print(tree)
# print()
# print()
# print()
# print(tree.strSequences(10, 15))

# tree = guide.runUPGMA(aln, 'fractional')
# guide.writeNewickFile('Files_WS2/my_cyp1a1_upgma.nwk', tree)

import csv as csv
import pandas as pd

def get_yeasts() -> list[str]:
  """ Get a list of all the yeast species """
  df = pd.read_csv("Files_WS2/sugars.csv")
  return list(df["Yeast"])

# def process_fa(path: str, whitelist: list[str]):
#   """
#   Process the fa file with yeast sequences.
#   We replace white space with _ for the label,
#   And exclude any species not in the list.
#   """
#   lines: list[str] = []

#   with open(path, 'rt') as fafile:
#     all_seq = list()
#     label = None
#     protein = ""
#     lines = fafile.readlines()

#   for line in lines:
#     line = line.strip()
#     if line.startswith(">"):
#       parts = line.split(" ")
#       name = parts[0][1:]
#       print(parts)

#       # not interested in sequences outside the yeast in sugars.csv.
#       if name not in whitelist:
#         continue

#       if label is not None:
#         data = {
#           "label": f"{label}",
#           "sequence": protein
#         }

#         # reset sequence for next one
#         protein = ""
#         all_seq.append(data)

#       label = f"{line.replace(" ", "_")} {parts[1]}"
#     else:
#       protein += line

#   all_seq.append({
#     "label": label,
#     "sequence": protein
#   })

#   return all_seq

# def clean_fasta():
#   """
#   Process the fa file and write to a new file, "select.fa"
#   """
#   yeasts = get_yeasts()

#   fa = process_fa("Files_WS2/MalS.fa", yeasts)
#   with open("select.fa", "a") as fafile:
#     for data in fa:
#       fafile.write(data["label"])
#       fafile.write("\n")
#       fafile.write(data["sequence"])
#       fafile.write("\n")


def exercise_7():
  """ 
  Curate the MalS.fa dataset
  """

  # Get all yeast from sugars.csv
  yeasts = get_yeasts()

  mals = guide.readFastaFile('Files_WS2/MalS.fa', guide.Protein_Alphabet)
  cnt = 0
  select = []

  # traverse mals dataset
  for seq in mals:

    # we are only interested in yeast also in the sugars.csv list
    if seq.name in yeasts:

      # white space -> _
      new_name = (seq.name + "_" + seq.annot).replace(" ", "_")
      seq.name = new_name
      cnt += 1
      select.append(seq)

  guide.writeFastaFile('select.fa', select)
  print('Selected', cnt, 'sequences')


# exercise_7()

# import part_b as part_b

def exercise_8():
  aln = guide.readClustalFile('exercise7.aln', guide.Protein_Alphabet)
  consensus = part_b.exercise_5(aln)
  # print("HERE", consensus)
  # index = [173, 231, 232, 233, 234, 294, 295, 324, 437]
  # print(consensus)
  # for i in index:
    # print(f"{i-1} -> {consensus[i-1]}")
    # print("Consensus Sequence is:\n\n", consensus)


ima = "".join(["MTISSAHPETEPKWWKEATFYQIYPASFKDSNDDGWGDMKGIASKLEYIKELGADAIWIS",
"PFYDSPQDDMGYDIANYEKVWPTYGTNEDCFALIEKTHKLGMKFITDLVINHCSSEHEWF",
"KESRSSKTNPKRDWFFWRPPKGYDAEGKPIPPNNWKSYFGGSAWTFDEKTQEFYLRLFCS",
"TQPDLNWENEDCRKAIYESAVGYWLDHGVDGFRIDVGSLYSKVVGLPDAPVVDKNSTWQS",
"SDPYTLNGPRIHEFHQEMNQFIRNRVKDGREIMTVGEMQHASDETKRLYTSASRHELSEL",
"FNFSHTDVGTSPLFRYNLVPFELKDWKIALAELFRYINGTDCWSTIYLENHDQPRSITRF",
"GDDSPKNRVISGKLLSVLLSALTGTLYVYQGQELGQINFKNWPVEKYEDVEIRNNYNAIK",
"EEHGENSEEMKKFLEAIALISRDHARTPMQWSREEPNAGFSGPSAKPWFYLNDSFREGIN",
"VEDEIKDPNSVLNFWKEALKFRKAHKDITVYGYDFEFIDLDNKKLFSFTKKYNNKTLFAA",
"LNFSSDATDFKIPNDDSSFKLEFGNYPKKEVDASSRTLKPWEGRIYISE"])

i = 158 - 1
print(ima[i-5:i+5])
print(ima)

mine = "".join(["---MTISSAHPETEPKWWKEATIYQIYPASFKDSN-----------NDGWGDLKGIASKLEYIKELGVDA",
"IWICPFYDSPQDDMGYDIANYEKVWPTYGTNEDCFALIEKTHKLGMKFITDLVINHCSSEHEWFKESRSS",
"KTNPKRDWFFWRPPKGYDAEGKPIPPNNWRSFFGGSAWTFDEKTQEFYLRLFASTQPDLNWENEDCRKAI",
"YESAVGYWLDHGVDGFRIDVGSLYSKVPGLPDAPVTDENSKWQHSDPFTMNGPRIHEFHQEMNKFMRNRV",
"-KDGREIMTVGEVQHGSDETKRLYTSASRHELSELFNFSHTDVGTSPKFRYNLVPFELKDWKVALAELFR",
"FINGTDCWSTIYLENHDQPRSITRFGDDSPKNRVISGKLLSVLLVSLTGTLYVYQGQELGQINF-KNWPI",
"EKYEDVEVRNNYKAIKEEHGENSK---EMKKFLEGIALISRDHARTPMPWTKEEPNAGFSG---PDAKPW",
"FYLNESFREGINAEDESKDPNSVLNFWKEALQFRKAHKDITVYGYDFEFIDLDNKKLFSFTKK--Y-DNK",
"TLFAALNFSSDEIDFTIPNDSASFKLEFGNYPDKEVDASSRTLKPWEGRIYISE--"])

# print("\nmine", mine)
# exercise_8()
# def exercise_9():
#   """ """

def exercise_10():
  aln = guide.readClustalFile('exercise7.aln', guide.Protein_Alphabet)
  phy = guide.runUPGMA(aln, 'poisson')
  # print(phy) # print so we can copy into jstree 
  return phy

# exercise_10()

def exercise_11():
  """ Annotate tree with metabolized sugars """
  df = pd.read_csv("Files_WS2/sugars.csv")
  yeasts = list(df["Yeast"])                        
  aln = guide.readClustalFile('exercise7.aln', guide.Protein_Alphabet)
  phy = guide.runUPGMA(aln, 'poisson')
  df = df.set_index(df['Yeast'])

  for col_name in df.columns[1:]:
    phy_copy = str(phy)
    t = 0
    f = 0
    for yeast in yeasts:
      can_metabolize = df[col_name][yeast]
      # add the labels using string replacement
      phy_copy = phy_copy.replace(yeast, f"{yeast} {can_metabolize}")
      # debugging - can see how many sugars can be metabolized
      if can_metabolize: 
        t += 1
      else:
        f += 1

    # Prints: <Melizitose> => 8 / 14 = 57.14285714285714
    # Useful for debugging
    if col_name == "Melizitose":
      print(phy_copy)
    print(f"{col_name} => {t} / {f + t} = {(t / (t + f)) * 100}")


# exercise_11()

aln = guide.readClustalFile('exercise7.aln', guide.Protein_wGAP)
# tree.putAlignment(aln)
# tree.parsimony()
# print(tree)
print()
print()
print()
# print(tree.strSequences(10, 15))

def getConsensusForColumn(aln, colidx) -> str:
    symcnt = {}
    for seq in aln.seqs:
        mysym = seq[colidx]
        if mysym in symcnt:
            symcnt[mysym] += 1
        else:
            symcnt[mysym] = 1
    consensus = ""
    maxcnt = 0
    for mysym in symcnt:
        if symcnt[mysym] > maxcnt:
            maxcnt = symcnt[mysym]
            consensus = mysym
    return consensus


def exercise_12():
    cols = len(aln[0])
    con_seq = ""
    # loop each col and grab the next char
    for i in range(cols):
        char = getConsensusForColumn(aln, i)
        # append it
        con_seq += char
    # return a Sequence instance
    return guide.Sequence(con_seq, guide.Protein_wGAP)

exercise_12()