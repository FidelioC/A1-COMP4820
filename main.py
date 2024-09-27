from Bio import SeqIO
from Bio.Seq import Seq
import click
import re


def read_fasta_file():
    seq_record = SeqIO.read("Sorangium_cellulosum_19lines.fasta", "fasta")
    simple_seq = Seq("ATGA")

    return seq_record


def brute_force(sequence_text, sequence_list_patterns):
    index = 0
    count = 0
    text_length = len(sequence_text)

    while index < text_length:
        # iterating the list of patterns
        for sequence_pattern in sequence_list_patterns:
            pattern_length = len(sequence_pattern)
            if (index + pattern_length < text_length) and (
                sequence_text.seq[index : (index + pattern_length)] == sequence_pattern
            ):
                count += 1
        index += 1
    return count


@click.command()
@click.option("--pattern", multiple=True)
def commands_processing(pattern):
    """command lines processing"""

    sequence_text = read_fasta_file()
    sequence_list_patterns = pattern
    print(sequence_list_patterns)
    print(brute_force(sequence_text, sequence_list_patterns))


if __name__ == "__main__":
    commands_processing()
