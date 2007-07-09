#!/usr/local/bin/perl

#File Names
$pwd = `pwd`;
chomp($pwd);
$output_filename = $pwd . "/StandardShaders.cc";

#Get Input File List
`rm -f $output_filename`;
@files = `find ./StandardShaders -maxdepth 2 -type f \! -path "*.svn*"`;

# Create Output File
open(OUT, ">$output_filename");

# Print Initial Code
print OUT "\#include <vw/GPU/StandardShaders.h>\n\n";
print OUT "namespace vw { namespace GPU {\n\n";
print OUT "std::map<std::string, char*> standard_shaders_map;\n\n";
print OUT "void init_standard_shaders() {\n";

#Iterate Through Input Files, printing map insert functions
foreach $filename (@files) {
  chomp($filename);
  open(IN, $filename);
  print OUT "standard_shaders_map[\"", $filename, "\"] = \" \\\n";
  while ($line = <IN>) {
    chomp($line);
    print OUT $line, " \\\n";
  }
  print OUT "\";\n\n";
}

# Print End Code
print OUT "}\n";
print OUT "} } // namespace GPU, namespace vw"
