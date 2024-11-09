#!/usr/bin/env python3


import subprocess
import sys

#____________________________ Step 1

def list_formation(blast_output:str) -> list[list[str]]:
    rows = blast_output.splitlines() # split text on \n
    list = []
    for row in rows:
        item = row.split("\t") # item is a list of items split by \t
        list.append(item)
    #print(list)
    return list


# getting full length hits (qlen==length) && pident > 80
def step_one(primer_file: str, assembly_file: str) -> list[list[str]]:
    cmd = " ".join([
    "blastn",
    "-query", primer_file,
    "-subject", assembly_file,
    "-task", "blastn-short",
    "-outfmt", "'6 std qlen'", "|", 
    "awk", "'$3>80 && $4==$13'", "|",
    "sort", "-n", "-k", "9,9"
    ])

    # capturing the output
    result = subprocess.run(cmd,capture_output=True, text=True, shell=True, executable="/bin/bash")
    # print(cmd)
    output = result.stdout
    return list_formation(output)



#____________________________ Step 2


def are_primers_facing2(start1, stop1, start2, stop2, max_amplicon):
    """Check if the primers are facing each other and within the max amplicon size."""
    if (start2 > stop2 and start1 < stop1): # if start is reverse orientation and stop in reverse orientation
        return (stop1 < stop2 and stop2-stop1 <= max_amplicon) 
    elif (start2 < stop2 and start1 > stop1): # if start in forward orientation and stop in reverse orientation
        return (stop2 < stop1 and stop1-stop2 <= max_amplicon)
    else:
        False



def step_two(sorted_hits: list[list[str]], max_amplicon_size: int) -> list[tuple[list[str]]]:
    primer_pairs = []
    seen = set()
    for hit in sorted_hits:
        for every_other_hit in sorted_hits:
            start_current = int(hit[8])
            stop_current = int(hit[9])
            start_other = int(every_other_hit[8])
            stop_other = int(every_other_hit[9])


            if are_primers_facing2(start_current, stop_current, start_other, stop_other, max_amplicon_size):
                pair = (hit, every_other_hit) 
                hit_tup = tuple(hit)
                every_other_tup = tuple(every_other_hit)
                # pair_tuple = (hit_tup, every_other_tup) # not hashable bc tried to use a list but theyre mutable so cant be an element in a set
                pair_tuple = tuple(sorted((hit_tup, every_other_tup), key=lambda x: x[0]))

                # pair_sorted = sorted(pair_tuple, key=lambda x: x[0]) # convert inner lists to tuples and sort so all pairs are in format (F, R)
                # sorted returns a list
                # pair_sorted = tuple(pair_sorted)
        
                if pair_tuple not in seen:  # check if the pair is unique
                    seen.add(pair_tuple)
                    primer_pairs.append(pair)
                else:
                    continue
            else:
                continue

        else:
            continue

    
    return primer_pairs



#____________________________ Step 3

def step_three(hit_pairs: list[tuple[list[str]]], assembly_file: str) -> str:
    start = ()
    # input is a list of tuples; inside each tuple is a pair of lists

    f = open("test2.bed", "w")
                
    for pair in hit_pairs: # inside list of tuples
        contig = pair[0][1]
        coords_1 = list((pair[0][8],pair[0][9]))
        coords_1 = [int(coord) for coord in coords_1]
        amplicon_start = max(coords_1)
        coords_2 = list((pair[1][8],pair[1][9]))
        coords_2 = [int(coord) for coord in coords_2]
        amplicon_end = min(coords_2) - 1
        #print(amplicon_start, amplicon_end)
        f.write(f"{contig}\t{amplicon_start}\t{amplicon_end}\n")

    f.close()

    cmd = " ".join([
    "seqtk",
    "subseq", assembly_file,
    f.name
    ])

    # capturing the output
    result = subprocess.run(cmd,capture_output=True, text=True, shell=True, executable="/bin/bash")
    output = result.stdout
    return output



#____________________________ Accumulated function run

def ispcr(primer_file: str, assembly_file: str, max_amplicon_size: int) -> str:
    sorted_list = step_one(primer_file, assembly_file)
    hit_pairs = step_two(sorted_list, max_amplicon_size)
    amplicons = step_three(hit_pairs, assembly_file)
    return amplicons


