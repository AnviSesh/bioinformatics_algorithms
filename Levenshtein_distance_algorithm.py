from Bio import SeqIO
from io import StringIO
import numpy as np


def read_sequence(prompt: str) -> str:
    mode = input(f"{prompt} — Enter 'F' for filename or 'P' to paste FASTA: ").strip().upper()

    if mode == "F":
        fname = input("Path to FASTA file: ").strip()
        record = next(SeqIO.parse(fname, "fasta"))
    elif mode == "P":
        print("Paste FASTA content. Type END on a new line when finished:")
        lines = []
        while True:
            line = input()
            if line.strip() == "END":
                break
            lines.append(line)
        fasta = "\n".join(lines)
        record = next(SeqIO.parse(StringIO(fasta), "fasta"))
    else:
        raise ValueError("Invalid option. Use 'F' or 'P'.")

    return str(record.seq)


def levenshtein_distance(s: str, t: str, show_matrix: bool = False) -> int:
    m, n = len(s), len(t)
    dp = np.zeros((m + 1, n + 1), dtype=int)

    for i in range(m + 1):
        dp[i, 0] = i
    for j in range(n + 1):
        dp[0, j] = j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if s[i - 1] == t[j - 1] else 1
            dp[i, j] = min(
                dp[i - 1, j] + 1,
                dp[i, j - 1] + 1,
                dp[i - 1, j - 1] + cost,
            )

    if show_matrix:
        print("\nDP matrix:")
        print(dp)

    return dp[m, n]


def main():
    seq1 = read_sequence("First sequence")
    seq2 = read_sequence("Second sequence")

    dist = levenshtein_distance(seq1, seq2, show_matrix=True)
    print("\nLevenshtein distance:", dist)


if __name__ == "__main__":
    main()
