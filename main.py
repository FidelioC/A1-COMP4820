from Bio import SeqIO
import click


# ========= BRUTE FORCE ========== #
def brute_force(sequence_text, sequence_list_patterns, dictionary):
    done_count = 0
    index = 0
    count = 0
    text_length = len(sequence_text)

    # going over the text one by one
    while done_count < len(sequence_list_patterns):
        # iterating the list of patterns
        for sequence_pattern in sequence_list_patterns:
            pattern_length = len(sequence_pattern)
            # check if they're the same pattern or if it's already went over the pattern length
            if index + pattern_length > text_length:
                done_count += 1
            elif (
                sequence_text.seq[index : (index + pattern_length)] == sequence_pattern
            ):
                count += 1
                dictionary[sequence_pattern].append(index)
        index += 1
    return dictionary


# ========= BMH ========== #
def shift_table_text(sequence_text, sequence_pattern):
    """create shift table for text (genome) pre-processing"""
    set_text = set(sequence_text)
    pattern_length = len(sequence_pattern)
    shift_table = {}
    index = pattern_length - 2  # ignore the last char

    # init text shift table with len(pattern)
    for item in set_text:
        shift_table[item] = pattern_length

    # update the shift table step by step, iterating the sequence_pattern
    while index >= 0:
        current_char = sequence_pattern[index]
        if (
            shift_table[current_char] == pattern_length
            and shift_table[current_char] != None
        ):
            shift_table[current_char] = pattern_length - index - 1

        index -= 1

    return shift_table


def shift_table_pattern(sequence_pattern):
    """create shift table for pattern pre-processing"""


def boyer_moore_horspool(sequence_text, sequence_list_patterns, dictionary):
    dict_indices = {}  # dictionary to safe the current index for each pattern
    done_count = 0  # used to see if all sequence patterns have reached the end
    index = 0
    count = 0
    text_length = len(sequence_text)

    # init all indices at 0 (at the start)
    for pattern in sequence_list_patterns:
        sequence_list_patterns[pattern] = 0

    # going over the text one by one
    while done_count < len(sequence_list_patterns):
        # iterating the list of patterns
        for sequence_pattern in sequence_list_patterns:
            pattern_length = len(sequence_pattern)
            # check if they're the same pattern or if it's already went over the pattern length
            if index + pattern_length > text_length:
                done_count += 1
            elif (
                sequence_text.seq[index : (index + pattern_length)] == sequence_pattern
            ):
                count += 1
                dictionary[sequence_pattern].append(index)
        index += 1
    return dictionary


# ========= Others ========== #
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
def commands_processing(pattern):
    # init
    sequence_text = read_fasta_file()
    sequence_list_patterns = pattern
    dictionary = initialize_output_dictionary(sequence_list_patterns)

    # algorithms
    dictionary = brute_force(sequence_text, sequence_list_patterns, dictionary)

    # output
    print_output(dictionary)


if __name__ == "__main__":
    commands_processing()

    # print(shift_table_text("auroraisauroraiguess", "aurora"))
