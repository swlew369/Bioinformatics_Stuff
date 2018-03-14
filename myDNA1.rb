#!/usr/bin/env ruby
# This Ruby script reads a DNA sequence and
# creates a codon count for all 3 Reading Frames of the DNA sequence
require 'bio'

# create a DNA sequence
seq = Bio::Sequence::NA.new("atggccattgaaggacat")
seq_1 = seq[1,seq.length]
seq_2 = seq[2,seq.length]
# translate to protein
prot = seq.translate

# prove that it worked
puts
puts seq   # => "atggccattgaatga"
puts seq_1
puts seq_2
puts prot  # => "MAIE*"
puts
# Generates a sample 100bp sequence.
seq1 = Bio::Sequence::NA.new("aatgacccgt" * 10)

# Naming this sequence as "testseq" and print in FASTA format
# (folded by 60 chars per line).
puts seq1.to_fasta("testseq", 50)
puts
codon_usage1 = Hash.new(0)
seq.window_search(3, 3) { |s| codon_usage1[s] += 1 }
puts "Codon Frequency for Reading Frame One:"
#puts codon_usage1
codon_usage1.each { |key, value| puts "#{key} #{value}" }
puts
codon_usage2 = Hash.new(0)
seq_1.window_search(3, 3) { |s| codon_usage2[s] += 1 }
puts "Codon Frequency for Reading Frame Two:"
#puts codon_usage2
codon_usage2.each { |key, value| puts "#{key} #{value}" }
puts
codon_usage3 = Hash.new(0)
seq_2.window_search(3, 3) { |s| codon_usage3[s] += 1 }
puts "Codon Frequency for Reading Frame Three:"
#puts codon_usage3
codon_usage3.each { |key, value| puts "#{key} #{value}" }
puts
