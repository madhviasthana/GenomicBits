from bit_encoding import *
from python_libraries import *
from data_loading import *
data_path="/Users/madhviasthana/Documents/Bioninformatics Git Repo/task-1-hamming-distance-histogram-madhviasthana/data/queries.txt"
import argparse
from hamminghist import *

# Define file paths
db_path = "/Users/madhviasthana/Documents/Bioninformatics Git Repo/task-1-hamming-distance-histogram-madhviasthana/data/database.txt"
query_path = "/Users/madhviasthana/Documents/Bioninformatics Git Repo/task-1-hamming-distance-histogram-madhviasthana/data/queries.txt"


# 2 Bit-Encoding
kmer_length =  get_kmer_length(data_path)
encoded_array = efficient_bit_encoding(data_path, kmer_length)
print(encoded_array.shape)

# Task 1: Hamming Distance
process_histogram_from_encoded(encoded_array, encoded_array, 20)
