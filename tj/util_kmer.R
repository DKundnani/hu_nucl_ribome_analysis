# utilities for working with kmers

kmer.NUCS <- c("A", "C", "G", "T")

kmer.NUCS_PLOT <- c("A", "G", "C", "T")

kmer.get_num_kmers <- function(kmer_size) {
  4 ^ kmer_size
}

kmer.ALL_KMERS <- list()

kmer.get_all_kmers <- function(kmer_size) {
  Biostrings::mkAllStrings(
    alphabet = kmer.NUCS,
    width = kmer_size
  )
}

kmer.reverse_complement <- function(x) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(x)))
}
