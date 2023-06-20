import sys
import random

# List of cell barcodes
BC = []

# read all cell barcodes from file and put into BC
with open("../data/samples/SI-TT-H4/outs/filtered_feature_bc_matrix/barcodes.tsv", "r") as bc_file:
    line = bc_file.readline()
    
    while line:
        BC.append(line.strip())
        line = bc_file.readline()

# randomly sample from BC
BC = random.sample(BC, k=900)

"""
Generate unique celltag strings.

counts: int; number of celltags to generate
version: str; should be v1, v2 or v3. Library version of celltag to generate.
"""
def generate_unique_strings(count, version):
    v = ("v1", "v2", "v3")
    start = ("GGT", "GTGATG", "TGTACG")[v.index(version)]
    unique_strings = set()
    while len(unique_strings) < count:
        middle_part = ''.join(random.choices('AGTC', k=8))
        unique_string = f'{start}{middle_part}GAATTC'
        unique_strings.add(unique_string)
    return list(unique_strings)

# create celltags for all libraries
v1 = generate_unique_strings(900, "v1")
v2 = generate_unique_strings(900, "v2")
v3 = generate_unique_strings(900, "v3")

# create weights for each celltag for each library
v1_w = [random.gauss(mu=30, sigma=20) for _ in range(900)]
v2_w = [random.gauss(mu=30, sigma=20) for _ in range(900)]
v3_w = [random.gauss(mu=30, sigma=20) for _ in range(900)]
bc_w = [random.gauss(mu=30, sigma=20) for _ in range(900)]

# sample barcodes and celltags based on weight for each library version
v1_BC = random.choices(BC, k=3500, weights=bc_w)
v1_seq = random.choices(v1, weights=v1_w, k=3500)
v2_BC = random.choices(BC, k=3000, weights=bc_w)
v2_seq = random.choices(v2, weights=v2_w, k=3000)
v3_BC = random.choices(BC, k=2000, weights=bc_w)
v3_seq = random.choices(v3, weights=v3_w, k=2000)


#v2_BC.extend(v3_BC)
#v1_BC.extend(v2_BC)
#v2_seq.extend(v3_seq)
#v1_seq.extend(v2_seq)

# add barcodes and celltags seperatly
BC = v1_BC + v2_BC +v3_BC
SEQ = v1_seq + v2_seq + v3_seq


# create needed UB information for each created sequence
UB = set()
while len(UB) < 8500:
    ub = ''.join(random.choices('AGTC', k=12))
    UB.add(ub)
UB = list(UB)

# zip the different columns together
sam_seq = list(zip(BC, SEQ, UB))


# add data to existing .sam file in .sam file format.
with open("../data/bam/data.sam", "a") as sam:
    for i, seq in enumerate(SEQ):
        bc = BC[i]#elem[0]
        cr = bc.split("-")[0]
        #seq = elem[1]
        ub = UB[i]#elem[2]
        
        q = "F"*len(seq)
        line = f"A00642:316:H7FHCDRX2:1:2177:27959:22451	4	*	0	0	*	*	0	0	{seq}	{q}	NH:i:0	HI:i:0	AS:i:37	nM:i:0	uT:A:1	ts:i:30	RG:Z:SI-TT-H4:0:1:H7FHCDRX2:1	xf:i:0	CR:Z:{cr}	CY:Z:FFFF:FFFFFFFFFFF	CB:Z:{bc}	UR:Z:{ub}	UY:Z:FFFFFFFFFFFF	UB:Z:{ub}\n"
        sam.write(line)
        