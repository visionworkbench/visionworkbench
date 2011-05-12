#!/usr/bin/env perl
# __BEGIN_LICENSE__
# Copyright (C) 2006-2011 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


use strict;
use warnings;
use File::Basename qw/dirname/;

#File Names
my $srcdir = dirname($0);
my $builddir = `pwd`;
chomp($builddir);

my $output_filename = "${builddir}/Shaders.cc";

#Get Input File List
`rm -f $output_filename`;
my @files = `find ${srcdir}/Shaders -maxdepth 2 -type f \! -path "*.svn*"`;

# Create Output File
open(OUT, '>', $output_filename);

# Print Initial Code
print OUT "\#include <vw/GPU/Shaders.h>\n\n";
print OUT "namespace vw { namespace GPU {\n\n";
print OUT "std::map<std::string, const char*> standard_shaders_map;\n\n";
print OUT "void init_standard_shaders() {\n";

#Iterate Through Input Files, printing map insert functions
foreach my $filename (@files) {
  chomp($filename);
  $filename =~ m#.*/Shaders/(.*)#;
  my $mapname = $1;
  open(IN, $filename);
  print OUT "standard_shaders_map[\"", $mapname, "\"] = \" \\\n";
  while (my $line = <IN>) {
    chomp($line);
    print OUT $line, " \\\n";
  }
  print OUT "\";\n\n";
}

# Print End Code
print OUT "}\n";
print OUT "} } // namespace GPU, namespace vw\n"
