# Genome Assembly with Overlap Graphs

## Overview

This project implements a *de novo* genome assembly pipeline using overlap graphs. *De novo* assembly aims to reconstruct a complete genome sequence from a set of fragmented reads without relying on a reference genome. This pipeline focuses on exploring the impact of various parameters on assembly performance, including read length, read count, minimum overlap, and error correction strategies.

The core of the assembly process is a greedy overlap-layout-consensus approach that constructs contigs (contiguous sequences) by iteratively merging reads based on their overlapping regions.

## Features

*   **Read Generation:** Simulates DNA reads from a reference genome (default: PhiX-174) with configurable read length, read count, and error rate.
*   **Overlap Graph Construction:** Builds an overlap graph to identify potential adjacencies between reads based on sequence overlap. Utilizes multiprocessing for parallel computation of overlaps.
*   **Greedy Assembly:** Implements a greedy algorithm to assemble reads into contigs by iteratively merging the best overlapping pairs.
*   **Error Correction:** Includes a k-mer frequency-based filtering method to correct errors in reads before assembly.
*   **Assembly Quality Assessment:** Calculates key assembly statistics such as N50, total assembled length, and correctness (alignment to the reference genome).
*   **Parameter Sensitivity Analysis:** Generates plots to visualize the impact of read length, read count, and error rates on assembly quality.

## Dependencies

*   Python 3.x
*   Biopython (`pip install biopython`)
*   Matplotlib (`pip install matplotlib`)
*   Tqdm (`pip install tqdm`)

## Installation

1.  Clone the repository:

    ```
    git clone 
    cd genomeAssembly
    ```

2.  Install the required Python packages:

    ```
    pip install biopython matplotlib tqdm
    ```

## Usage

1.  **Prepare the reference genome:** Ensure that the `NC_001422.fasta` file (PhiX-174 genome) is located in a directory named `files` in the same directory as the script.
    ```
    mkdir files
    # Download the reference genome (example using wget)
    wget -O files/NC_001422.fasta "https://ftp.ncbi.nlm.nih.gov/genomes/all/ ভাইরাস / бактериофаг_phiX_174_sensu_lato_uid13879/NC_001422.fasta"
    ```

2.  **Run the main script:**

    ```
    python TestPhixGraphs.py
    ```


## Running Tests and Creating Graphs

After the program has finished running, it saves all the requested plots as `.png` images.

## Parameters

The following parameters can be adjusted in the `main()` function of the script:

*   `N`: Number of reads.
*   `l`: Read length.
*   `min_overlap`: Minimum overlap length required to merge reads.
*   `error_rate`: Simulated error rate in the generated reads.
*   `k`: k-mer length for error correction.
*   `frequency_threshold`: Minimum frequency threshold for k-mer filtering.

## K-mer Based Filtering

The assembly pipeline utilizes k-mer based filtering for error correction. Here's a breakdown:

1.  **K-mer Generation:** All k-mers of length 'k' are generated from each read.

2.  **K-mer Counting:** The frequency of each k-mer across all reads is counted.

3.  **Filtering:** Reads are filtered based on whether all their k-mers occur at least `frequency_threshold` times. Reads that contain infrequent k-mers (likely containing errors) are discarded.

Careful tuning of `k` and `frequency_threshold` is required to balance error correction effectiveness with the risk of removing genuine reads.

## Output

The script prints assembly statistics to the console, including:

*   Number of contigs
*   Longest contig length
*   Total assembled length
*   N50
*   Correctness (percentage of aligned bases to the reference genome)

The script also generates plots visualizing the relationship between assembly parameters and contig length, saved as `.png` images:

*   `longest_contig_vs_N_error_free.png`: Longest contig length vs. number of reads (N) in error-free reads.
*   `longest_contig_vs_L_error_free.png`: Longest contig length vs. read length (L) in error-free reads.
*   `Total_contig_vs_N_with_and_without_ec.png`: Total contig length vs. number of reads (N) for different error rates, with and without error correction (k-mer filtering).

## Contributing

Please feel free to submit pull requests with bug fixes, improvements, or new features.

Special thanks to Sigal Tsabari for genetic consulting.
## License

Copyright (c) 2025 Yotam Tsabari

