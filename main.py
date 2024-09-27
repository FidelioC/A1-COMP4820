from Bio import SeqIO


def read_fasta_file():
    for seq_record in SeqIO.parse("Sorangium_cellulosum.fasta", "fasta"):
        print(seq_record.id)
        print(repr(seq_record.seq))
        print(len(seq_record))


if __name__ == "__main__":
    read_fasta_file()
