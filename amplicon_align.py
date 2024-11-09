#!/usr/bin/env python3
import argparse
import magnumopus

'''
Usage: ./amplicon_align.py -1 data/Pseudomonas_aeruginosa_PAO1.fna -2 data/Pseudomonas_protegens_CHA0.fna -p data/rpoD.fna -m 2000 --match 1 --mismatch=-1 --gap=-1
'''

def parse_arguments():
    parser = argparse.ArgumentParser(description='Align amplicons from two assembly files.')
    parser.add_argument('-1', '--assembly1', required=True, help='Path to the first assembly file')
    parser.add_argument('-2', '--assembly2', required=True, help='Path to the second assembly file')
    parser.add_argument('-p', '--primer', required=True, help='Path to the primer file')
    parser.add_argument('-m', '--max_amplicon_size', type=int, required=True, help='Maximum amplicon size')
    parser.add_argument('--match', type=int, default=1, help='Score for sequence match (default: 1)')
    parser.add_argument('--mismatch', type=int, default=-1, help='Score for sequence mismatch (default: -1)')
    parser.add_argument('--gap', type=int, default=-1, help='Penalty for gap (default: -1)')
    
    return parser.parse_args()


def reverse_complement(string):
    complement_dict = {"A": "T", "T":"A", "G":"C", "C":"G"}
    reverse_complement = ""
    for letter in string:
        if letter in complement_dict.keys():
            reverse_complement += complement_dict[letter]
        else:
            reverse_complement += letter
    reverse_complement_final = reverse_complement[::-1] # use splicing to reverse the string
    return reverse_complement_final



def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Perform isPCR to find amplicons in both assembly files
    amplicon1 = magnumopus.ispcr(args.primer, args.assembly1, args.max_amplicon_size)
    amplicon2 = magnumopus.ispcr(args.primer, args.assembly2, args.max_amplicon_size)

    amplicon1 = amplicon1.split('\n', 1)[1].strip()
    amplicon2 = amplicon2.split('\n', 1)[1].strip()

    
    if amplicon1 is None or amplicon2 is None:
        print("No valid amplicons found in one or both assembly files.")
        return

    # Align the amplicons using Needleman-Wunsch algorithm
    # Assuming that there is only one amplicon from each file so don't have to worry about multiple
    aln1, score1 = magnumopus.needleman_wunsch(amplicon1, amplicon2, args.match, args.mismatch, args.gap)

    # Reverse complement the second amplicon for possible better alternative alignment
    amplicon2_rc = reverse_complement(amplicon2)
    aln2, score2 = magnumopus.needleman_wunsch(amplicon1, amplicon2_rc, args.match, args.mismatch, args.gap)

    # Determine the best alignment (higher score)
    if score1 >= score2:
        best_aln = aln1
        best_score = score1
    else:
        best_aln = aln2
        best_score = score2

    # Print the best alignment and its score
    print("Best Alignment:\n")
    print("\n".join(best_aln))
    print("\nAlignment Score:", best_score)

if __name__ == '__main__':
    main()
