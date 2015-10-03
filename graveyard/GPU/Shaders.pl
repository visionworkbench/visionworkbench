#!/usr/bin/env perl
# __BEGIN_LICENSE__
#  Copyright (c) 2006-2013, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NASA Vision Workbench is licensed under the Apache License,
#  Version 2.0 (the "License"); you may not use this file except in
#  compliance with the License. You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
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
