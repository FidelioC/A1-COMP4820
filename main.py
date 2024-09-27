from Bio import SeqIO
from Bio.Seq import Seq
import re


def read_fasta_file():
    seq_record = SeqIO.read("Sorangium_cellulosum_19lines.fasta", "fasta")
    simple_seq = Seq("ATGA")

    return seq_record


def brute_force(sequence_text, sequence_pattern):
    index = 0
    count = 0
    pattern_length = len(sequence_pattern)
    text_length = len(sequence_text)

    while (index + pattern_length) < text_length:
        # print(sequence_text.seq[index : (index + pattern_length)])
        if sequence_text.seq[index : (index + pattern_length)] == sequence_pattern:
            count += 1
        index += 1
    return count


if __name__ == "__main__":
    sequence_pattern = "ATGA"
    sequence_text = read_fasta_file()
    print(brute_force(sequence_text, sequence_pattern))
