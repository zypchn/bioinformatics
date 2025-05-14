import numpy as np

seq1 = "TGTTACGG"  
seq2 = "GGTTGACTA"   

n1 = len(seq1) + 1    # rows
n2 = len(seq2) + 1    # cols

matrix = np.zeros((n1, n2))

gap_penalty = -2
match = +1
mismatch = -1

def init_mx():
    for i in range(1, n1):
        for j in range(1, n2):
            ismatch = seq1[i-1] == seq2[j-1]
            match_score = match if ismatch else mismatch
            matrix[i][j] = max((matrix[i-1][j] + gap_penalty), (matrix[i][j-1] + gap_penalty), (matrix[i-1][j-1] + match_score), 0)
    return matrix
        
def get_max_positions(matrix):
    max_list = (np.where(matrix == matrix.max()))
    max_idx = []
    for i in range(len(max_list[0])):
        idx = max_list[0][i], max_list[1][i]
        max_idx.append(idx)
    return max_idx

def print_stars(seq1, seq2):
    stars = ""
    for i in range(len(seq1)):
        if (seq1[i] == seq2[i]):
            stars += "*"
        else:
            stars += " "
    print(stars)

matrix = init_mx()
max_pos = get_max_positions(matrix)

sample_pos = max_pos[0]
row = sample_pos[0]
col = sample_pos[1]

align_1 = []
align_2 = []
align_1.append(seq1[row-1])
align_2.append(seq2[col-1])
row -= 1
col -=1

traverse = True
while traverse:
    score = matrix[row][col]
    
    if (score > 0):
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
                
    traverse = score > 0             


align_seq1 = "".join(align_1)[::-1]
align_seq2 = "".join(align_2)[::-1]
print(align_seq1)
print_stars(align_seq1, align_seq2)
print(align_seq2)