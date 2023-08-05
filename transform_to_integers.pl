#!/usr/bin/perl
use strict;
use warnings;

# Open the input file,MUST be in the same folder as the script
open(my $input_fh, "<", "countdataFloat") or die "Can't open input file: $!";

# Open the output file
open(my $output_fh, ">", "countdata") or die "Can't open output file: $!";


my $header1 = <$input_fh>;
my $header2 = <$input_fh>;
print $output_fh $header1;
print $output_fh $header2;

# Read through the input file line by line
while (my $line = <$input_fh>) {
  chomp($line);

  # Split the line into fields
  my @fields = split(/\t/, $line);

  # Print the first field (gene ID) unchanged
  print $output_fh "$fields[0]\t";

  # Convert the remaining fields to integers and print them
  for (my $i = 1; $i < scalar(@fields); $i++) {
    my $value = int($fields[$i] + 0.5);
    print $output_fh "$value\t";
  }

  print $output_fh "\n";
}

# Close the input and output files
close($input_fh);
close($output_fh);
