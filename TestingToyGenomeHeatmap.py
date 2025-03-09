import multiprocessing as mp
from Bio import SeqIO
import random
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from collections import Counter,defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

mp.set_start_method('fork', force=True)

# Function to compute the maximum overlap between two reads
def max_overlap(read1, read2, min_overlap=1):
    max_possible = min(len(read1), len(read2))
    curr_overlap = 0
    for i in range(min_overlap, max_possible + 1):
        if read1[-i:] == read2[:i]:
            curr_overlap = i
    return curr_overlap


# Function to compute overlaps for a single read (used in parallel)
def compute_overlap(read_index, reads, min_overlap):
    overlaps = {}
    read = reads[read_index]
    for j, second_read in enumerate(reads):
        if read_index != j:
            curr_overlap = max_overlap(read, second_read, min_overlap)
            if curr_overlap:
                overlaps[j] = curr_overlap
    return read_index, overlaps


# Function to build the suffix array (same as before)
def build_suffix_array(s):
    return sorted(range(len(s)), key=lambda i: s[i:])


# Function to build the LCP array (same as before)
def build_lcp_array(s, sa):
    n = len(s)
    rank = [0] * n
    for i in range(n):
        rank[sa[i]] = i
    lcp = [0] * n
    k = 0
    for i in range(n - 1):
        if rank[i] == 0:
            continue
        j = sa[rank[i] - 1]
        while i + k < n and j + k < n and s[i + k] == s[j + k]:
            k += 1
        lcp[rank[i]] = k
        if k > 0:
            k -= 1
    return lcp


# Function to find overlaps using the suffix and LCP arrays (same as before)
def find_overlaps(sa, lcp, read_lengths, min_overlap):
    overlaps = defaultdict(dict)
    n = len(sa)
    for i in range(1, n):
        if lcp[i] >= min_overlap:
            read1, pos1 = next((j, sa[i - 1] - sum(read_lengths[:j])) for j in range(len(read_lengths)) if
                               sa[i - 1] < sum(read_lengths[:j + 1]))
            read2, pos2 = next((j, sa[i] - sum(read_lengths[:j])) for j in range(len(read_lengths)) if
                               sa[i] < sum(read_lengths[:j + 1]))
            if read1 != read2:
                if pos1 == 0 and pos2 + lcp[i] == read_lengths[read2]:
                    overlaps[read2][read1] = lcp[i]
                elif pos2 == 0 and pos1 + lcp[i] == read_lengths[read1]:
                    overlaps[read1][read2] = lcp[i]
    return overlaps


# Function to build the overlap graph in parallel
def build_overlap_graph_parallel(reads, min_overlap=1):
    overlap_graph = {}
    with mp.Pool() as pool:
        results = pool.starmap(compute_overlap, [(i, reads, min_overlap) for i in range(len(reads))])

        for read_index, overlaps in results:
            overlap_graph[read_index] = overlaps

    return overlap_graph

# Function for the greedy assembly algorithm



def find_best_overlap_parallel(graph, min_overlap):
    optimal_overlap = 0
    best_pair = (None, None)
    for read1 in graph:
        for read2, curr_overlap in graph[read1].items():
            if curr_overlap > optimal_overlap:
                best_pair = (read1, read2)
                optimal_overlap = curr_overlap
    return optimal_overlap, best_pair


def greedy_assembly_overlap_graph(reads, min_overlap=1):
    reads = reads.copy()
    total_iterations = len(reads)

    with ProcessPoolExecutor() as executor:
        with tqdm(total=total_iterations, desc="Assembling reads") as pbar:
            while len(reads) > 1:
                graph = build_overlap_graph_parallel(reads, min_overlap)

                # Parallelize the search for the best overlap
                chunk_size = max(1, len(graph) // executor._max_workers)
                futures = []
                for i in range(0, len(graph), chunk_size):
                    chunk = dict(list(graph.items())[i:i + chunk_size])
                    futures.append(executor.submit(find_best_overlap_parallel, chunk, min_overlap))

                optimal_overlap = 0
                best_pair = (None, None)
                for future in as_completed(futures):
                    curr_overlap, curr_pair = future.result()
                    if curr_overlap > optimal_overlap:
                        optimal_overlap = curr_overlap
                        best_pair = curr_pair
                pbar.update(1)
                pbar.set_postfix({"Reads remaining": len(reads)})
                if optimal_overlap < min_overlap:
                    break

                read1, read2 = best_pair
                new_read = reads[read1] + reads[read2][optimal_overlap:]
                reads = [read for i, read in enumerate(reads) if i not in (read1, read2)]
                reads.append(new_read)



    return reads


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

def run_assembly_and_get_total_contig_length(genome, N, l, error_rate, kmer_filter=False, k=10, frequency_threshold=2, min_overlap=20):
    """Runs the assembly pipeline and returns the total length of the contigs."""

    # Generate reads
    reads = generate_reads(genome, N, l, error_rate)

    # K-mer filtering
    if kmer_filter:
        filtered_reads = kmer_based_filtering(reads, k, frequency_threshold)
    else:
        filtered_reads = reads

    # Assemble filtered reads
    contigs = greedy_assembly_overlap_graph(filtered_reads, min_overlap)

    # Calculate total contig length
    total_contig_length = sum(len(c) for c in contigs)

    return total_contig_length

def create_toy_genome(length=1000):
    """Creates a random toy genome of specified length."""
    return ''.join(random.choice('ATCG') for _ in range(length))

def generate_heatmap_data(genome_length = 1000):
    """Generates the data for the heatmap."""

    N_values = [100, 500, 1000, 1500, 2000]
    L_values = [50, 100, 150, 200, 250]
    error_rates = [0.00, 0.01, 0.02, 0.04, 0.08]
    min_overlaps = [5, 15, 25, 35, 45]

    results = {}

    # Use the new function to make the genome
    toy_genome = create_toy_genome(genome_length)

    for error_rate in error_rates:
        for min_overlap in min_overlaps:
            heatmap_data = np.zeros((len(N_values), len(L_values)))
            for i, N in enumerate(N_values):
                for j, L in enumerate(L_values):
                    total_contig_length = run_assembly_and_get_total_contig_length(toy_genome, N, L, error_rate, min_overlap=min_overlap)
                    heatmap_data[i, j] = total_contig_length
            results[(error_rate, min_overlap)] = heatmap_data

    return results, N_values, L_values, error_rates, min_overlaps

def plot_heatmaps(results, N_values, L_values, error_rates, min_overlaps):
    """Plots the heatmaps."""
    num_rows = len(error_rates)
    num_cols = len(min_overlaps)

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_cols * 5, num_rows * 4))

    for i, error_rate in enumerate(error_rates):
        for j, min_overlap in enumerate(min_overlaps):
            ax = axes[i, j] if num_rows > 1 and num_cols > 1 else (axes[i] if num_cols == 1 else axes[j])
            data = results[(error_rate, min_overlap)]
            im = ax.imshow(data, cmap="YlOrRd", origin="lower", aspect="auto")

            # Set ticks and labels
            ax.set_xticks(np.arange(len(L_values)))
            ax.set_yticks(np.arange(len(N_values)))
            ax.set_xticklabels(L_values)
            ax.set_yticklabels(N_values)

            # Rotate the tick labels and set their alignment.
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

            ax.set_title(f"p={error_rate:.2f}, min_overlap={min_overlap}")
            if j == 0:
                ax.set_ylabel("N")
            if i == len(error_rates) - 1:
                ax.set_xlabel("l")

    fig.tight_layout()
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.75)  # Add a colorbar
    plt.savefig("assembly_heatmap_toy_genome.png")
    plt.close()

# Main program
def main():
    # Set main genome size for toy genome
    genome_length = 1000

    # Load PhiX genome
    phix_genome = next(SeqIO.parse("files/NC_001422.fasta", "fasta")).seq

    # Parameters
    N = 400  # Number of reads
    l = 150  # Read length
    min_overlap = 20  # Minimum overlap for assembly
    error_rate = 0.02 #introduce error rate

    #Make small toy genome
    toy_genome = create_toy_genome(genome_length)

    # Generate reads
    reads = generate_reads(toy_genome, N, l,error_rate)
    # Apply k-mer based filtering
    k = 10  # k-mer length
    frequency_threshold = 2  # Minimum frequency to keep a read
    filtered_reads = kmer_based_filtering(reads, k, frequency_threshold)

    print(f"Original number of reads: {len(reads)}")
    print(f"Number of reads after filtering: {len(filtered_reads)}")

    # Calculate average coverage (using filtered reads)
    genome_length = len(toy_genome)
    avg_coverage = (len(filtered_reads) * l) / genome_length
    print("genome length:", genome_length)
    print(f"Average coverage: {avg_coverage:.2f}x")

    # Assemble filtered reads
    contigs = greedy_assembly_overlap_graph(filtered_reads, min_overlap)

    # Calculate assembly statistics
    N50 = calculate_N50(contigs)
    correctness = calculate_correctness(contigs, toy_genome)

    # Print results
    print("\nError rate:", error_rate)
    print(f"Number of contigs: {len(contigs)}")
    print(f"Longest contig length: {max(len(c) for c in contigs)}")
    print(f"Total assembled length: {sum(len(c) for c in contigs)}")
    print(f"N50: {N50}")
    print(f"Correctness: {correctness:.2f}%")

    # Generate and plot the heatmap
    results, N_values, L_values, error_rates, min_overlaps = generate_heatmap_data()
    plot_heatmaps(results, N_values, L_values, error_rates, min_overlaps)
    print("Heatmap created successfully!")

if __name__ == "__main__":
    main()
