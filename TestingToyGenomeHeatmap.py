import multiprocessing as mp
from Bio import SeqIO
import random
import matplotlib.pyplot as plt
from tqdm import tqdm
from collections import Counter

mp.set_start_method('fork', force=True)


def generate_toy_genome(size=1000):
    bases = ['A', 'T', 'C', 'G']
    genome = ''.join(random.choices(bases, k=size))

    return genome


def max_overlap(read1, read2, min_overlap=1, error_rate=0.02):
    """
    Calculates the maximum overlap between two reads, allowing for mismatches,
    but prioritizing perfect matches.
    """
    max_possible = min(len(read1), len(read2))
    mismatches = 0
    for i in range(max_possible, min_overlap - 1, -1):
        if read1[-i:] == read2[:i]:
            return i  # Perfect overlap found

    for i in range(max_possible, min_overlap - 1, -1):
        if i < max_possible:
            if read1[-i - 1] != read2[i - 1]:
                mismatches -= 1
            if read1[-i] != read2[0]:
                mismatches += 1
        else:
            mismatches = sum(c1 != c2 for c1, c2 in zip(read1[-i:], read2[:i]))

        if mismatches / i <= error_rate:
            return i
    return 0  # No acceptable overlap found


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
    """
    Calculates the correctness of the assembly
    """
    ref_len = len(reference_genome)
    coverage = [0] * ref_len

    for contig in contigs:
        contig_len = len(contig)
        if contig_len > ref_len:
            continue


        ref_hash = hash(reference_genome[:contig_len])
        contig_hash = hash(contig)

        for i in range(ref_len - contig_len + 1):
            if contig_hash == ref_hash:

                if contig == reference_genome[i:i + contig_len]:
                    for j in range(contig_len):
                        coverage[i + j] = 1

            if i < ref_len - contig_len:
                ref_hash = hash(reference_genome[i + 1:i + contig_len + 1])

    total_aligned = sum(coverage)
    correctness = (total_aligned / ref_len) * 100
    return correctness

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

def run_assembly_and_get_total_contig_length(N, l, error_rate, kmer_filter=False, k=10, frequency_threshold=2,min_overlap=20):
    """Runs the assembly pipeline and returns the total length of the contigs."""
    toy_genome = generate_toy_genome()

    # Generate reads
    reads = generate_reads(toy_genome, N, l, error_rate)

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

def create_hash_map():
    """Creates and populates the hash map with assembly results."""
    toy_genome = generate_toy_genome()
    error_rates = [0.00, 0.01, 0.03, 0.05]
    min_overlaps = [5, 15, 25, 40]
    read_lengths = [50, 100, 150, 200]
    num_reads_values = [100, 200, 400,600]

    results_map = {}

    for error_rate in error_rates:
        results_map[error_rate] = {}
        for min_overlap in min_overlaps:
            results_map[error_rate][min_overlap] = {}
            for read_length in read_lengths:
                results_map[error_rate][min_overlap][read_length] = {}
                for num_reads in num_reads_values:
                    total_contig_length = run_assembly_and_get_total_contig_length(
                        num_reads, read_length, error_rate, min_overlap=min_overlap
                    )
                    results_map[error_rate][min_overlap][read_length][num_reads] = total_contig_length

    return results_map

def create_combined_heatmap(results_map):
    """Generates a single image containing all heatmaps in a grid."""
    error_rates = [0.0, 0.01, 0.03, 0.05]
    min_overlaps = [5, 15, 25, 40]
    read_lengths = [50, 100, 150,200]
    num_reads_values = [100, 200, 400,600]

    # Determine the number of rows and columns for the subplots
    num_rows = len(error_rates)
    num_cols = len(min_overlaps)

    # Create a single figure with all subplots
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(20, 20))

    # Loop through the error rates and min overlaps to create the heatmaps
    for i, error_rate in enumerate(error_rates):
        for j, min_overlap in enumerate(min_overlaps):
            # Prepare data for the heatmap
            heatmap_data = []
            for num_reads in num_reads_values:
                row = [results_map[error_rate][min_overlap][read_length][num_reads] for read_length in read_lengths]
                heatmap_data.append(row)

            # Select the appropriate subplot
            ax = axes[i, j]

            # Create the heatmap
            im = ax.imshow(heatmap_data, cmap="YlOrRd", origin="lower")

            # Set labels
            ax.set_xticks(np.arange(len(read_lengths)))
            ax.set_yticks(np.arange(len(num_reads_values)))
            ax.set_xticklabels(read_lengths)
            ax.set_yticklabels(num_reads_values)

            # Add title
            ax.set_title(f"p={error_rate}, min_overlap={min_overlap}")
            ax.set_xlabel("l")
            ax.set_ylabel("N")

            # Rotate the tick labels and set their alignment.
            plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                     rotation_mode="anchor")

            # Ensure the subplots are properly spaced
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout to make room

            # Create colorbar
            if i == 0 and j == num_cols - 1:  # Add colorbar to the last subplot
                cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.95)

    # Add a title to the entire figure
    fig.suptitle("Total Contig Length Heatmaps", fontsize=16)

    # Save the combined figure
    plt.savefig("combined_heatmaps.png")
    plt.close()

# Main program
def main():
    
    toy_genome = generate_toy_genome()
    # Parameters
    N = 200  # Number of reads
    l = 150  # Read length
    min_overlap = 20  # Minimum overlap for assembly
    error_rate = 0.02 #introduce error rate

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

    # Create and populate the hash map
    results_map = create_hash_map()

    # Create the graphs
    create_combined_heatmap(results_map)

if __name__ == "__main__":
    main()
