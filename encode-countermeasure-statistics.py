import math
import itertools
from random import randbytes
from ascii_graph import Pyasciigraph

def count_ones(bytearray):
    count = 0
    for byte in bytearray:
        binary_string = bin(byte)[2:].zfill(8)  # Convert byte to binary string with leading zeros
        count += binary_string.count('1')
    return count

def random_bits_distribution(n):
    ret = []
    for i in range(n):
        r = count_ones(randbytes(32))
        ret.append(r)
    return ret

def print_histogram(data):
    graph = Pyasciigraph()
    for line in graph.graph(None, data):
        print(line)

dist = random_bits_distribution(500000)

histodata_bits = []
histodata_permutations = []
for i in sorted(set(dist)):
    histodata_bits.append((f"{i}", dist.count(i)))
    histodata_permutations.append((math.floor(math.log2(math.comb(256, i))), dist.count(i)))

print_histogram(histodata_permutations)

