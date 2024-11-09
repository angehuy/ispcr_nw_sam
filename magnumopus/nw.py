#!/usr/bin/env python3

import subprocess
import sys

def make_matrix(x_dim: int, y_dim: int) -> list[list]:
    matrix = []
    for i in range(x_dim):
        matrix.append(y_dim * [0])
    return matrix

def populate_initial(matrix: list[list], seq_a_len: int, seq_b_len: int, seq_a, seq_b):
    count_1 = 0
    count_2 = -1

    for i, row in enumerate(matrix): # enumerate is to keep track of index while iterating over iterable
       
        if i == 0:
            for j in range(seq_b_len): # has to be seq_b_len for it to work
                matrix[i][j] = count_1
                count_1 -= 1

        else: # when not in the first
            for j in range(seq_b_len):
                if j == 0:
                    matrix[i][0] = count_2
                    count_2 -= 1
    return matrix
        

def calculate_score(matrix: list[list], seq_a_len: int, seq_b_len: int, seq_a:str, seq_b:str, match: int, mismatch: int, gap: int):
    for i in range(1, seq_a_len): # skip the first row
        for j in range(1, seq_b_len): # skip the first column
            score_diagonal = matrix[i-1][j-1] + (match if seq_a[i-1] == seq_b[j-1] else mismatch)
            score_left = matrix[i][j-1] + gap
            score_up = matrix[i-1][j] + gap
            matrix[i][j] = max(score_diagonal, score_left, score_up)
    # print(matrix)
    return matrix

# Diagonal move means store both sequences
# Vertical move means store gap in the horizontal sequence
# Horizontal move means store gap in the vertical sequence
def traceback(matrix, seq_a_len, seq_b_len,seq_a, seq_b, match, mismatch, gap):
    i,j = seq_a_len-1, seq_b_len-1 # set i and j to the lengths so we can traceback at the end
    aligned_a, aligned_b = '', ''
    score = 0
    while i>0 and j>0:
        if i >0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] + (match if seq_a[i-1] == seq_b[j-1] else mismatch):
            aligned_a += seq_a[i-1]
            aligned_b += seq_b[j-1]
            i-= 1
            j -=1
        elif i > 0 and matrix[i][j] == matrix[i-1][j] + gap: # veritcal move
            aligned_a += seq_a[i-1] # store in horizontal sequence
            aligned_b += "-"
            i -= 1
        else:
            aligned_a += "-"
            aligned_b += seq_b[j-1]
            j -= 1
    aligned_a_final = aligned_a[::-1]
    aligned_b_final = aligned_b[::-1]

    score = 0
    for a, b in zip(aligned_a_final, aligned_b_final):
        if a == b:
            score += match
        elif a == '-' or b == '-':
            score += gap
        else:
            score += mismatch

    return (aligned_a_final, aligned_b_final), score




def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    seq_a_len = len(seq_a) + 1 # add 1 for the gap row
    seq_b_len = len(seq_b) + 1 # add 1 for the gap column
    matrix = make_matrix(seq_a_len, seq_b_len)
    init_matrix = populate_initial(matrix, seq_a_len, seq_b_len, seq_a, seq_b)
    calculated_matrix = calculate_score(init_matrix, seq_a_len, seq_b_len, seq_a, seq_b, match, mismatch, gap)
    aligned_sequences, score = traceback(calculated_matrix, seq_a_len, seq_b_len, seq_a, seq_b, match, mismatch, gap)
    return aligned_sequences, score

