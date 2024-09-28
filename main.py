from Bio import SeqIO
import click
import bmh
import brute_force
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ========= main ========== #
def read_fasta_file():
    seq_record = SeqIO.read("Sorangium_cellulosum_19lines.fasta", "fasta")

    return seq_record


def initialize_output_dictionary(sequence_list_patterns):
    """
    Init output directory with empty list
    Purpose: store 'count' and 'location' of the patterns
    """
    dictionary = {}
    for pattern in sequence_list_patterns:
        if pattern not in dictionary:
            dictionary[pattern] = []
    return dictionary


def print_output(dictionary):
    # print(dictionary)
    for item in dictionary:
        print(
            f"Pattern: {item} \n\tCount: {len(dictionary[item])} \n\tLocations: {dictionary[item]}"
        )


@click.command()
@click.option("--pattern", multiple=True)
@click.option("--algorithm")
def commands_processing(pattern, algorithm):
    # init
    sequence_text = read_fasta_file()
    sequence_list_patterns = list(pattern)
    dictionary = initialize_output_dictionary(sequence_list_patterns)

    # algorithms
    if algorithm == "bmh":
        print("BMH: ")
        all_tables = bmh.build_all_shift_tables(sequence_text, sequence_list_patterns)
        dictionary = bmh.boyer_moore_horspool(
            sequence_text, sequence_list_patterns, dictionary, all_tables
        )
    else:
        print("BRUTE FORCE: ")
        dictionary = brute_force.brute_force(
            sequence_text, sequence_list_patterns, dictionary
        )

    # output
    print_output(dictionary)


if __name__ == "__main__":
    commands_processing()
