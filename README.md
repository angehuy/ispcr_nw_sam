**magnumopus module** contains packages **ispcr** and **nw**, performing in silico PCR (isPCR), a simulation of the
polymerase chain reaction (PCR) commonly used in wet-lab settings to amplify target DNA sequences, and the Needleman-Wunsch algorithm, respectively.

">" denotes a `.py` file intended to be run in the terminal.

**ispcr.py** performs 3 steps:
1) Perform a BLASTN search of a query file (primers) against the subject file (assembly) to identify primer annealing locations.
2) Identify pairs of primer annealing sites that would yield a valid amplicon (primers facing each other and within close proximity).
3) Extract amplified sequences within the assembly using `seqtk`.

> `q1.py` takes as input the required files for `ispcr.py` functions.  
--> **Usage:** `./q1.py`

**nw.py** performs sequence alignment based on the Needleman-Wunsch algorithm:  
1) Creates an alignment matrix.  
2) Populates the matrix with gap values.  
3) Calculates matrix scores based on sequence matches and mismatches.  
4) Performs traceback to determine potential alignments.

> `q2.py` takes as input 2 random DNA sequences.  
--> **Usage:** `./q2.py`

**amplicon_align.py** accepts command-line arguments.  
--> **Usage:**  
`./amplicon_align.py -1 data/Pseudomonas_aeruginosa_PAO1.fna -2 data/Pseudomonas_protegens_CHA0.fna -p data/rpoD.fna -m 2000 --match 1 --mismatch=-1 --gap=-1`

**Steps performed by amplicon_align.py**:  
1) Performs isPCR to find amplicons in both assembly files (assuming 1 amplicon per assembly file).  
2) Aligns the amplicons using the Needleman-Wunsch algorithm.  
3) Reverse complements the second amplicon to explore possible better alternative alignments.  
4) Determines the best alignment (indicated by a higher score) and prints the alignment along with its score.
