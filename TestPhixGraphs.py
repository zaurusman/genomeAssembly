import multiprocessing as mp
from Bio import SeqIO
import random
import matplotlib.pyplot as plt
from tqdm import tqdm

mp.set_start_method('fork', force=True)

def max_overlap(read1, read2, min_overlap=1):
    max_possible = min(len(read1), len(read2))
    # Search from largest to smallest
    for i in range(max_possible, min_overlap-1, -1):
        if read1.endswith(read2[:i]):
            return i
    return 0

def compute_overlap(read_index, reads, min_overlap):

    overlaps = {}
    read = reads[read_index]
    for j, second_read in enumerate(reads):
        if read_index != j:
            curr_overlap = max_overlap(read, second_read, min_overlap)
            if curr_overlap:
                overlaps[j] = curr_overlap
    return read_index, overlaps

def build_overlap_graph_parallel(reads, min_overlap=1, num_processes=None):
    if num_processes is None:
        num_processes = mp.cpu_count()

    overlap_graph = {}
    with mp.Pool(num_processes) as pool:
        results = pool.starmap(compute_overlap, [(i, reads, min_overlap) for i in range(len(reads))])

        for read_index, overlaps in results:
            overlap_graph[read_index] = overlaps

    return overlap_graph

def introduce_errors(read, error_rate):

    mutated_read = list(read)
    for base_index in range(len(read)):
        if random.random() < error_rate:

            possible_bases = ['A', 'T', 'C', 'G']
            current_base = read[base_index]
            possible_bases.remove(current_base)
            mutated_read[base_index] = random.choice(possible_bases)
    return "".join(mutated_read)

def generate_reads(genome, num_reads, read_length, error_rate=0.0):

    genome_length = len(genome)
    reads = []
    for _ in range(num_reads):
        start = random.randint(0, genome_length - read_length)
        read = genome[start:start + read_length]
        if error_rate > 0:
            read = introduce_errors(read, error_rate)
        reads.append(str(read))
    return reads

def greedy_assembly_overlap_graph(reads, min_overlap=1):
    reads = reads.copy()
    total_iterations = len(reads)

    with tqdm(total=total_iterations, desc="Assembling reads") as pbar:
        while True:
            optimal_overlap = 0
            best_pair = (None,None)
            graph = build_overlap_graph_parallel(reads, min_overlap)

            for read1 in graph:
                for read2, curr_overlap in graph[read1].items():
                    if curr_overlap > optimal_overlap:
                        best_pair = [read1, read2]
                        optimal_overlap = curr_overlap
            pbar.update(1)
            pbar.set_postfix({"Reads remaining": len(reads)})
            if optimal_overlap < min_overlap:
                break

            new_reads = [reads[i] for i in range(len(reads)) if i not in best_pair]
            new_read = reads[best_pair[0]] + reads[best_pair[1]][optimal_overlap:]
            new_reads.append(new_read)
            reads = new_reads


    return reads

def calculate_N50(contigs):

    contig_lengths = sorted([len(c) for c in contigs], reverse=True)
    total_length = sum(contig_lengths)
    cumulative_length = 0
    for length in contig_lengths:
        cumulative_length += length
        if cumulative_length >= total_length / 2:
            return length
    return 0  # Should not happen if contigs is not empty

def calculate_correctness(contigs, reference_genome):
    """Calculates the correctness of the assembly by aligning contigs to the reference genome (assuming we have one for this assignment)"""
    total_aligned_length = 0
    for contig in contigs:
        # Find the best alignment of the contig to the reference genome
        best_alignment_length = 0
        for i in range(len(reference_genome) - len(contig) + 1):
            alignment_length = 0
            for j in range(len(contig)):
                if contig[j] == reference_genome[i + j]:
                    alignment_length += 1
            best_alignment_length = max(best_alignment_length, alignment_length)
        total_aligned_length += best_alignment_length
    return total_aligned_length / len(reference_genome)*100

def kmer_based_filtering(reads, k, frequency_threshold):
    """Filter reads based on k-mer frequencies."""
    kmer_counts = Counter()
    for read in reads:
        kmers = generate_kmers(read, k)
        kmer_counts.update(kmers)
    filtered_reads = []
    for read in reads:
        kmers = generate_kmers(read, k)
        if all(kmer_counts[kmer] >= frequency_threshold for kmer in kmers):
            filtered_reads.append(read)
    return filtered_reads

def generate_kmers(read, k):
    """Generate k-mers from a read."""
    return [read[i:i+k] for i in range(len(read) - k + 1)]

from collections import Counter

def run_assembly_and_get_total_contig_length(N, l, error_rate, kmer_filter=False, k=10, frequency_threshold=2,min_overlap=20):
    """Runs the assembly pipeline and returns the total length of the contigs."""
    # Load PhiX genome
    phix_genome = next(SeqIO.parse("files/NC_001422.fasta", "fasta")).seq

    # Generate reads
    reads = generate_reads(phix_genome, N, l, error_rate)

    # K-mer filtering
    if kmer_filter:
        filtered_reads = kmer_based_filtering(reads, k, frequency_threshold)
    else:
        filtered_reads = reads

    # Assemble filtered reads
    contigs = greedy_assembly_overlap_graph(filtered_reads, min_overlap)
    #if there are no contigs
    if not contigs:
        return 0

    # Calculate total contig length
    total_contig_length = sum(len(c) for c in contigs)

    return total_contig_length

def create_graphs(min_overlap):
    """Creates the graphs for assembly quality assessment."""

    # Load PhiX genome
    phix_genome = next(SeqIO.parse("files/NC_001422.fasta", "fasta")).seq

    # 1. Longest contig length with changing N (error-free)
    N_values = [100, 200, 300, 400, 600, 800]
    longest_contig_lengths_N = []
    for N in N_values:
        total_contig_length = run_assembly_and_get_total_contig_length(N, 150, 0.0)
        longest_contig_lengths_N.append(total_contig_length)

    plt.figure(figsize=(10, 6))
    plt.plot(N_values, longest_contig_lengths_N, marker='o', linestyle='-', color='blue')
    plt.title("Longest Contig Length vs. N (Error-Free)")
    plt.xlabel("Number of Reads (N)")
    plt.ylabel("Longest Contig Length")
    plt.grid(True)
    plt.savefig("longest_contig_vs_N_error_free.png")
    plt.close()

    # 2. Longest contig length with changing L (error-free)
    L_values = [50, 100, 150, 200, 250,300]
    longest_contig_lengths_L = []
    for l in L_values:
        total_contig_length = run_assembly_and_get_total_contig_length(400, l, 0.0,min_overlap=min_overlap)
        longest_contig_lengths_L.append(total_contig_length)

    plt.figure(figsize=(10, 6))
    plt.plot(L_values, longest_contig_lengths_L, marker='o', linestyle='-', color='green')
    plt.title("Longest Contig Length vs. L (Error-Free)")
    plt.xlabel("Read Length (L)")
    plt.ylabel("Longest Contig Length")
    plt.grid(True)
    plt.savefig("longest_contig_vs_L_error_free.png")
    plt.close()

    # 3. Longest contig length with changing N for different error values (with and without error correction)
    error_rates = [0.0, 0.01, 0.02, 0.04]
    colors = ['red', 'orange', 'purple', 'green']  # Matching colors for each error value

    plt.figure(figsize=(12, 8))

    for i, error_rate in enumerate(error_rates):
        # Without error correction
        longest_contig_lengths_no_ec = []
        for N in N_values:
            total_contig_length = run_assembly_and_get_total_contig_length(N, 150, error_rate, kmer_filter=False,min_overlap=min_overlap)
            longest_contig_lengths_no_ec.append(total_contig_length)
        plt.plot(N_values, longest_contig_lengths_no_ec, marker='o', linestyle='-', color=colors[i], label=f"Error Rate {error_rate} (No EC)")

        # With k-mer based error correction
        longest_contig_lengths_ec = []
        for N in N_values:
            total_contig_length = run_assembly_and_get_total_contig_length(N, 150, error_rate, kmer_filter=True,min_overlap=min_overlap)
            longest_contig_lengths_ec.append(total_contig_length)
        plt.plot(N_values, longest_contig_lengths_ec, marker='x', linestyle='--', color=colors[i], label=f"Error Rate {error_rate} (With EC)")

    plt.title("Total Contig Length vs. N for Different Error Rates (With and Without Error Correction)")
    plt.xlabel("Number of Reads (N)")
    plt.ylabel("Total Contig Length")
    plt.grid(True)
    plt.legend()
    plt.savefig("Total_contig_vs_N_with_and_without_ec.png")
    plt.close()

    print("Graphs created successfully!")

# Main program
def main():
    # Load PhiX genome
    phix_genome = next(SeqIO.parse("files/NC_001422.fasta", "fasta")).seq

    # Parameters
    N = 400  # Number of reads
    l = 100  # Read length
    min_overlap = 20  # Minimum overlap for assembly
    error_rate = 0.00 #introduce error rate

    # Generate reads
    reads = generate_reads(phix_genome, N, l,error_rate)
    # Apply k-mer based filtering
    k = 6  # k-mer length
    frequency_threshold = 3  # Minimum frequency to keep a read
    filtered_reads = kmer_based_filtering(reads, k, frequency_threshold)

    print(f"Original number of reads: {len(reads)}")
    print(f"Number of reads after filtering: {len(filtered_reads)}")

    # Calculate average coverage (using filtered reads)
    genome_length = len(phix_genome)
    avg_coverage = (len(filtered_reads) * l) / genome_length
    print("genome length:", genome_length)
    print(f"Average coverage: {avg_coverage:.2f}x")

    # Assemble filtered reads
    contigs = greedy_assembly_overlap_graph(filtered_reads, min_overlap)

    # Calculate assembly statistics
    N50 = calculate_N50(contigs)
    correctness = calculate_correctness(contigs, phix_genome)

    # Print results
    print("\nError rate:", error_rate)
    print(f"Number of contigs: {len(contigs)}")
    print(f"Longest contig length: {max(len(c) for c in contigs)}")
    print(f"Total assembled length: {sum(len(c) for c in contigs)}")
    print(f"N50: {N50}")
    print(f"Correctness: {correctness:.2f}%")

    #create_graphs(min_overlap)  # Call the function to create the graphs

if __name__ == "__main__":
    main()
