from SWalign import init_score_mx as sw_score_mx
from SWalign import get_score as get_sw_score
from NWalign import init_score_mx as nw_score_mx
from NWalign import get_score as get_nw_score

match_score = 1
mismatch_score = -1
gap_penalty = -2

seq1 = "TGTTACGG"
seq2 = "GGTTGACTA"

sw_matrix = sw_score_mx(seq1, seq2, match_score, mismatch_score, gap_penalty)
sw_score = get_sw_score(sw_matrix)


nw_matrix = nw_score_mx(seq1, seq2, match_score, mismatch_score, gap_penalty)
nw_score = get_nw_score(nw_matrix)

print(f"Sequence 1 : {seq1}")
print(f"Sequence 2 : {seq2}")
print(f"SW Score : {sw_score}")
print(f"NW Score : {nw_score}")