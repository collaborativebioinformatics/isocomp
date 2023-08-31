# std lib
import logging
# external dependencies
from Bio import SeqIO

logger = logging.getLogger(__name__)

__all__ = ['fasta_to_fastq']

def fasta_to_fastq(fasta_file: str) -> None:
    """Convert a fasta file to a fastq file. This function simply iterates 
    over a fasta file and prints a fastq file line for each line of the fasta 
    file

    Args:
        fasta_file (str): path to a fasta file
    """

    logger.debug(fasta_file)

    for seq in SeqIO.parse(fasta_file, "fasta"):
        seq.letter_annotations["solexa_quality"] = [40] * len(seq)
        print(seq.format("fastq"), end='')
