# ========= BMH ========== #
def shift_table(sequence, pattern):
    """create shift table for text (genome) or pattern pre-processing
    depending on the input"""
    set_text = set(sequence)
    pattern_length = len(pattern)
    shift_table = {}
    index = pattern_length - 2  # ignore the last char

    # init text shift table with len(pattern)
    for item in set_text:
        shift_table[item] = pattern_length

    # update the shift table step by step, iterating the sequence_pattern
    while index >= 0:
        current_char = pattern[index]
        if (
            shift_table[current_char] == pattern_length
            and shift_table[current_char] != None
        ):
            shift_table[current_char] = pattern_length - index - 1

        index -= 1

    return shift_table


def build_all_shift_tables_text(sequence_text, sequence_list_patterns):
    """create all shift tables for each pattern for text processing"""
    all_shift_tables = {}
    for pattern in sequence_list_patterns:
        # ignore duplicate pattern
        if pattern not in all_shift_tables:
            all_shift_tables[pattern] = shift_table(sequence_text, pattern)

    return all_shift_tables


def build_all_shift_tables_pattern(sequence_list_patterns):
    """create all shift tables for each pattern for pattern processing"""
    all_shift_tables = {}
    for pattern in sequence_list_patterns:
        # ignore duplicate pattern
        if pattern not in all_shift_tables:
            all_shift_tables[pattern] = shift_table(pattern, pattern)

    return all_shift_tables


def boyer_moore_horspool(
    sequence_text, sequence_list_patterns, dictionary, all_shift_tables
):
    dict_indices = {}  # dictionary to safe the current index for each pattern
    index = 0  # start at index 0
    text_length = len(sequence_text)

    # init all indices at 0 (at the start)
    for sequence_pattern in sequence_list_patterns:
        dict_indices[sequence_pattern] = 0

    # going over the text one by one
    while sequence_list_patterns:
        # iterating the list of patterns
        for sequence_pattern in sequence_list_patterns:
            index = dict_indices[sequence_pattern]
            pattern_length = len(sequence_pattern)
            # check if they're the same pattern or if it's already went over the pattern length
            if index + pattern_length <= text_length:
                if (
                    sequence_text.seq[index : (index + pattern_length)]
                    == sequence_pattern
                ):
                    # save the found location in the dictionary
                    dictionary[sequence_pattern].append(index)

                # shift the current pattern index according to the last compared char from the text
                next_char_shift = sequence_text.seq[(index + pattern_length - 1)]

                # conditional check if symbol in T is in the shift table or no
                if next_char_shift in all_shift_tables[sequence_pattern]:
                    # shift according to table
                    dict_indices[sequence_pattern] += all_shift_tables[
                        sequence_pattern
                    ][next_char_shift]
                else:
                    # decide what value to shift, here I used len(pattern) for shifting
                    dict_indices[sequence_pattern] += pattern_length

            else:
                sequence_list_patterns.remove(sequence_pattern)

        # print(f"dict_indices {dict_indices}")
    return dictionary
