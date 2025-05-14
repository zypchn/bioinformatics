import numpy as np

def init_score_mx(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-2):
    n1 = len(seq1) + 1
    n2 = len(seq2) + 1
    matrix = np.zeros((n1, n2))
    for i in range(1, n1):
        for j in range(1, n2):
            ismatch = seq1[i-1] == seq2[j-1]
            match = match_score if ismatch else mismatch_score
            matrix[i][j] = max((matrix[i-1][j] + gap_penalty), (matrix[i][j-1] + gap_penalty), (matrix[i-1][j-1] + match), 0)
    return matrix

def get_max_positions(matrix):
    max_list = (np.where(matrix == matrix.max())) if matrix.max() > 1 else []
    if max_list == []:
        return max_list
    max_idx = []
    for i in range(len(max_list[0])):
        idx = max_list[0][i], max_list[1][i]
        max_idx.append(idx)
    return max_idx


def align_sequences(matrix, positions, seq1, seq2):
    if positions == []:
        print("No aligmnents were found!")
        return
    
    alignments = []
    for sample_pos in positions:
        row = sample_pos[0]
        col = sample_pos[1]

        align_1 = []
        align_2 = []
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
        alignments.append(align)
    return alignments

def print_stars(seq1, seq2):
    stars = ""
    for i in range(len(seq1)):
        if (seq1[i] == seq2[i]):
            stars += "*"
        else:
            stars += " "
    print(stars)
    
def print_result(alignments):
    for alignment in alignments:
        align_1 = alignment[0]
        align_2 = alignment[1]
        align_seq1 = "".join(align_1)[::-1]
        align_seq2 = "".join(align_2)[::-1]
        print(align_seq2)
        print_stars(align_seq1, align_seq2)
        print(align_seq1)
        print("\n")
    
def main(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-2):
    score_mx = init_score_mx(seq1, seq2)
    max_positions = get_max_positions(score_mx)
    alignments = align_sequences(score_mx, max_positions, seq1, seq2)
    print_result(alignments)
    

if __name__ == "__main__":
    
    mismatch_score = -1
    match_score = 1
    gap_penalty = -2
    
    seq1 = "TGTTACGG"  
    seq2 = "GGTTGACTA" 
     
    main(seq1, seq2)
