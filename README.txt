Introduction

ReadMill is a command-line program for processing BAM format
(http://samtools.sourceforge.net) short-read DNA sequence data,
particularly unmapped reads.

BAM format is being used as a compact alternative to Fastq format for
storing unmapped reads. ReadMill comprises a set of established Fastq
tools ported to operate on BAM data.

The goal is to provide a simple command-line program that

1) performs reporting tasks such as making quality plots and counts

2) performs data manipiulation tasks such as filtering, partitioning,
   clipping and annotating reads prior to mapping

Data manipulation will be read-group aware and will add suitable
metadata to SAM headers to record what changes have been made.
