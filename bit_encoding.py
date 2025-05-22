from codebase.python_libraries import *
""" 
- Read DNA k-mers from a file
- Convert them into compact 2-bit encoded integers
"""

def get_kmer_length(file_path):
    with open(file_path, 'r') as f:
        first_line = f.readline()
        kmer = first_line.split()[0]
        return len(kmer)

 # Method1: 2 Bit Encoding of Kmers using ASCII values and masking them

def bit_encoding_using_loop(seq):
    """ seq: sequence from file"""
    encoded = 0
    for char in seq:
        base_bits = (ord(char) >> 1) & 0b11
        encoded = (encoded << 2) | base_bits 

    return encoded

# create a numpy array of all the sequences
def create_np_array(kmer,kmer_length):
    encoded_data = []
    # Choose dtype
    if kmer_length <= 8:
        dtype = np.uint16
    elif kmer_length <= 16:
        dtype = np.uint32
    else:
        dtype = np.uint64

    np_encoded_array = np.array(encoded_data, dtype=dtype)

    #print("NumPy array shape:", np_encoded_array.shape)

    return np_encoded_array

"""_summary_:
    The Disadvantage of using method 1 is that: For smaller dataset its okay, but if dataset is larger then going through loops is not the best thing one can do, its expensive and time taking
    """

# Method 2: Use Numpy for 2 Bit Encoding

def efficient_bit_encoding(file_path, kmer_length):
    """
    Encode DNA sequences from a file into 2-bit representations using ASCII bits 2 and 3.
    Parameters:
    kmer_length: int, length of each k-mer
    
    Returns:
    np.ndarray: Array of encoded values
    """
    sequences = np.genfromtxt(file_path, dtype=str, usecols=0)
    
    # Convert sequences to a NumPy array of ASCII values
    seq_bytes = np.array([np.frombuffer(seq.encode('ascii'), dtype=np.uint8) for seq in sequences])
    
    encoded_bases = (seq_bytes >> 2) & 0b11

    # Choose appropriate dtype based on kmer_length
    if kmer_length <= 8:
        dtype = np.uint16
    elif kmer_length <= 16:
        dtype = np.uint32
    else:
        dtype = np.uint64
    
    shifts = 2 * (kmer_length - 1 - np.arange(kmer_length, dtype=dtype))
    shifts = shifts.reshape(1, kmer_length)
    
    encoded_bases = encoded_bases.astype(dtype=dtype)
    
    shifted_bases = np.left_shift(encoded_bases, shifts)
    encoded = np.sum(shifted_bases, axis=1)
    return encoded.astype(dtype)