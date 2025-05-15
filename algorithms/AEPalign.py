import numpy as np
from tkinter import *
from SWalign import init_score_mx
from NWalign import align_sequences

def compute_norm_score(score, length, length_scal):
    return score / (length + length_scal)

def compute_iteration(x, norm_score):
    return x - (2 * norm_score)

#seq1 = "TGTTACGG"
#seq2 = "GGTTGACTA"
seq1 = "TGTT"
seq2 = "GGTT"

mx = init_score_mx(seq1, seq2)
vectorized_func = np.vectorize(compute_iteration)
score = mx.max()
norm_score = compute_norm_score(0.18, 6, 10)
res = compute_iteration(mx, norm_score)

score = compute_norm_score(3, 6, 10)
print(mx - 2 * score)