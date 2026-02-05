#!/usr/bin/env python3

from Bio import SeqIO


def read_fasta_sequence(fasta_file: str) -> str:
    """Read the first sequence from a FASTA file."""
    record = next(SeqIO.parse(fasta_file, "fasta"))
    return str(record.seq)


def build_kmer_index(sequence: str, k: int) -> dict:
    """
    Build an index mapping each k-mer to a sorted list
    of start positions (0-based).
    """
    index = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i : i + k]
        index.setdefault(kmer, []).append(i)
    return index


def binary_search_left(arr: list, target: int) -> int:
    """Find leftmost index of target using binary search."""
    left, right = 0, len(arr)
    while left < right:
        mid = (left + right) // 2
        if arr[mid] < target:
            left = mid + 1
        else:
            right = mid
    return left


def binary_search_right(arr: list, target: int) -> int:
    """Find rightmost index of target using binary search."""
    left, right = 0, len(arr)
    while left < right:
        mid = (left + right) // 2
        if arr[mid] <= target:
            left = mid + 1
        else:
            right = mid
    return left


def get_kmer_positions(index: dict, kmer: str) -> list:
    """
    Return all positions where the k-mer occurs
    using binary search boundaries.
    """
    if kmer not in index:
        return []

    positions = index[kmer]

    left = binary_search_left(positions, positions[0])
    right = binary_search_right(positions, positions[-1])

    return positions[left:right]


def main():
    fasta_path = input("Enter path to FASTA file: ").strip()
    k = int(input("Enter k (3 for codons): ").strip())

    sequence = read_fasta_sequence(fasta_path)
    index = build_kmer_index(sequence, k)

    print("Sequence length:", len(sequence))
    print("Unique k-mers:", len(index))

    kmer = input("Enter codon/k-mer to search: ").strip().upper()

    positions = get_kmer_positions(index, kmer)

    if not positions:
        print(f"{kmer} not found in sequence.")
    else:
        print(f"{kmer} occurs {len(positions)} times.")
        print("Positions (0-based):", positions)


if __name__ == "__main__":
    main()
