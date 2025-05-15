import numpy as np
from tkinter import *

"""
Needleman-Wunsch Implementation
"""

# score matrix computation
def init_score_mx(seq1, seq2, match_score, mismatch_score, gap_penalty):
    n1 = len(seq1) + 1
    n2 = len(seq2) + 1
    matrix = np.zeros((n1, n2))
    matrix[0] = [gap_penalty * i for i in range(n2)]
    matrix[:, 0] = [gap_penalty * i for i in range(n1)]
    for i in range(1, n1):
        for j in range(1, n2):
            ismatch = seq1[i-1] == seq2[j-1]
            match = match_score if ismatch else mismatch_score
            matrix[i][j] = max((matrix[i-1][j] + gap_penalty), (matrix[i][j-1] + gap_penalty), (matrix[i-1][j-1] + match))
    return matrix

def get_score(matrix):
    return matrix[-1][-1]

# traceback
def align_sequences(matrix, position, seq1, seq2):
    if position == []:
        print("No aligmnents were found!")
        return None
    
    row, col = position[0], position[1]
    align_1, align_2 = [], []
    align_1.append(seq1[row-1])
    align_2.append(seq2[col-1])
    row -= 1
    col -=1
        
    while matrix[row][col] > 0:

        nuc1 = seq1[row-1]
        nuc2 = seq2[col-1]
        
        if (nuc1 == nuc2):
            align_1.append(nuc1)
            align_2.append(nuc2)
            row -= 1
            col -= 1
            
        else:
            up = matrix[row-1][col]
            left = matrix[row][col-1]
            
            if left > up:   # add up's nuc
                align_1.append(nuc1)
                align_2.append("_")
                col -= 1
                
            else:
                align_1.append("_")
                align_2.append(nuc2)
                row -= 1
            
    align = [align_1, align_2]
    return align

# for readability of the results
def print_stars(seq1, seq2):
    stars = ""
    for i in range(len(seq1)):
        if (seq1[i] == seq2[i]):
            stars += "*"
        else:
            stars += " "
    return (stars)

# format results    
def get_result(alignment):
    formatted_results = []
    align_1, align_2 = alignment[0], alignment[1]
    align_seq1 = "".join(align_1)[::-1]
    align_seq2 = "".join(align_2)[::-1]
    result = f"""
    {align_seq2}\n
    {print_stars(align_seq1, align_seq2)}
    {align_seq1}
    """
    formatted_results.append(result)
    return formatted_results

# GUI display 
def display_results(results, score):
    master = Tk()
    master.title("Needleman-Wunsch Alignment")
    master.minsize(350, 100)
    for i in range(len(results)):
        Label(master, text=results[i]).grid(row=0, column=i)
    Label(master, text=f"Score = {score}").grid(column=len(results))
    mainloop()

# main function    
def nwalign(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-2):
    score_mx = init_score_mx(seq1, seq2, match_score, mismatch_score, gap_penalty)
    #alignments = align_sequences(score_mx, max_positions, seq1, seq2)
    #formatted_results = get_result(alignments)
    #display_results(formatted_results, score)


# runner
if __name__ == "__main__":
    
    mismatch_score = -1
    match_score = 1
    gap_penalty = -2
    
    seq1 = "AATCG"  
    seq2 = "AACG" 
     
    #nwalign(seq1, seq2)
    
    matrix = init_score_mx(seq1, seq2, match_score, mismatch_score, gap_penalty)
    #pos = get_max_positions(matrix)
    #print(pos)