class Node:
    def __init__(self) -> None:
        self.output = []
        self.goto = {}
        self.fail = None


class Trie:
    def __init__(self, sequence_list_patterns) -> None:
        self.root = Node()
        self.build_trie(sequence_list_patterns)
        print(self.root.goto)

    def build_trie(self, sequence_list_patterns):
        """build the trie from the paterns"""
        # iterate each pattern
        for pattern in sequence_list_patterns:
            # for each pattern, always start at the root node
            curr_node = self.root
            # iterate each char for that pattern
            for char in pattern:
                # if there's no goto arc for that node, create a new one
                if char not in curr_node.goto:
                    curr_node.goto[char] = Node()
                # go to the next child node
                curr_node = curr_node.goto[char]
            # finished creating nodes for one pattern, this is the output node
            curr_node.output.append(char)

    def build_fail_arcs():
        """build fail arcs"""

    def aho_corasick():
        """set pattern matching using aho corasick"""
