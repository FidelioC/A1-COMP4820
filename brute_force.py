from Bio import SeqIO


# ========= BRUTE FORCE ========== #
def brute_force(sequence_text, sequence_list_patterns, dictionary):
    index = 0
    text_length = len(sequence_text)

    # going over the text one by one
    while sequence_list_patterns:
        # iterating the list of patterns
        for sequence_pattern in sequence_list_patterns:
            pattern_length = len(sequence_pattern)
            # check if they're the same pattern or if it's already went over the pattern length
            if index + pattern_length <= text_length:
                if (
                    sequence_text.seq[index : (index + pattern_length)]
                    == sequence_pattern
                ):
                    dictionary[sequence_pattern].append(index)
            else:
                sequence_list_patterns.remove(sequence_pattern)
        index += 1
    return dictionary
