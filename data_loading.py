# This function loads the data from the file

def load_data(data_path):
    """Load k-mers from a text file and return a list of k-mers.
     
    Args:
        data_path: Path to the text file containing k-mers and counts.
    
    Returns:
        kmers: List of k-mer strings.
        kmer_length: Length of the k-mers (assuming all are equal length).
    """
    kmers = []
    kmer_length = None
    
    with open(data_path, "r") as file:
        for line in file:
            if not line.strip():
                continue 
            kmer, count = line.strip().split()
            kmers.append(kmer)
            # Set kmer_length from the first k-mer (assuming all are equal length)
            if kmer_length is None:
                kmer_length = len(kmer)
            # Optional: Validate that all k-mers have the same length
            elif len(kmer) != kmer_length:
                raise ValueError(f"K-mer {kmer} has length {len(kmer)}, expected {kmer_length}")
    
    return kmers, kmer_length


    
