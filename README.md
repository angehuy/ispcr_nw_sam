**magnumopus module** contains packages **ispcr** and **nw**, performing  in silico PCR (isPCR), a simulation of the
polymerase chain reaction (PCR) that is commonly used in wet-lab settings to amplify target DNA sequences, and the needleman-wunsch algorithm respectively.
">" denotes a .py file intended to be run in the terminal

**ispcr.py** performs 3 steps
1) blastn of a query file (primers) against the subject file (assembly) to identify primer annealing locations
2) identify pairs of primer annealing sites that would yield a valid amplicon (primers facing each other and within close proximity)
3) extracting amplified sequences within assembly with seqtk

> q1.py takes as input for ispcr.py functions required files --> usage: ./q1.py

**nw.py** performs sequence alignment based on the needleman-wunsch algorithm
1) creates matrix
2) populates matrix with gap values
3) calculate matrix scores based on sequence matches, mismatches
4) traceback to determine potential alignments

> q2.py takes as input 2 random DNA sequence --> usage: ./q2.py

> **amplicon_align.py** accepts command-line arguments --> usage: ./amplicon_align.py -1 data/Pseudomonas_aeruginosa_PAO1.fna -2 data/Pseudomonas_protegens_CHA0.fna -p data/rpoD.fna -m 2000 --match 1 --mismatch=-1 --gap=-1
1) performs isPCR to find amplicons in both assembly files (assuming 1 amplicon per assembly file)
2) align the amplicons using needleman-wunsch algorithm
3) reverse complement the second amplicon for possible better alternative alignment
4) determine the best alignment (indicated by a higher score) and prints alignment with score
