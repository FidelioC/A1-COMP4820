from Bio import SeqIO
import click


def read_fasta_file():
    seq_record = SeqIO.read("Sorangium_cellulosum_19lines.fasta", "fasta")

    return seq_record


def brute_force(sequence_text, sequence_list_patterns, dictionary):
    index = 0
    count = 0
    text_length = len(sequence_text)

    # going over char by char
    while index < text_length:
        # iterating the list of patterns
        for sequence_pattern in sequence_list_patterns:
            pattern_length = len(sequence_pattern)
            # check if they're the same pattern
            if (index + pattern_length < text_length) and (
                sequence_text.seq[index : (index + pattern_length)] == sequence_pattern
            ):
                count += 1
                dictionary[sequence_pattern].append(index)
        index += 1
    return dictionary


def initialize_dictionary(sequence_list_patterns):
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
def commands_processing(pattern):
    # init
    sequence_text = read_fasta_file()
    sequence_list_patterns = pattern
    dictionary = initialize_dictionary(sequence_list_patterns)

    # algorithms
    dictionary = brute_force(sequence_text, sequence_list_patterns, dictionary)

    # output
    print_output(dictionary)


if __name__ == "__main__":
    commands_processing()
