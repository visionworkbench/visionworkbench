#!/usr/bin/perl
use strict;
use warnings;

use File::Slurp;
use File::Basename;

my %comment = (
    ".ac"   => "dnl",
    ".am"   => "#",
    ".cc"   => "//",
    ".cg"   => "//",
    ".glsl" => "//",
    ".h"    => "//",
    ".i"    => "//",
    ".m4"   => "dnl",
    ".mak"  => "#",
    ".pl"   => "#",
    ".py"   => "#",
    ".tcc"  => "//",
);

# Read the license text from __DATA__
my @license = read_file(\*DATA);

# process each line given on stdin
foreach my $filename (<>) {
    chomp $filename;
    print "Processing $filename\n";

    # get the extension, and skip it if we don't know about it
    my (undef, undef, $ext) = fileparse($filename, qr/\.[^.]*/);
    next unless exists $comment{$ext};

    my $file = read_file($filename);

    # Remove a license header if it exists
    $file =~ s/^[^\n]*__BEGIN_LICENSE__.*__END_LICENSE__[^\n]*$//ms;

    # Remove all blank files from the top of the file
    while ($file =~ s/^\s*\n//) {};

    # prepend the license text, prepending the comment string to each line.
    # Also, separate the license header from content by two blank lines
    $file = $comment{$ext} . join($comment{$ext}, @license) . "\n\n" . $file;

    write_file($filename, $file);
}

__DATA__
 __BEGIN_LICENSE__

 Copyright (C) 2006 United States Government as represented by the
 Administrator of the National Aeronautics and Space Administration
 (NASA).  All Rights Reserved.

 Copyright 2006 Carnegie Mellon University. All rights reserved.

 This software is distributed under the NASA Open Source Agreement
 (NOSA), version 1.3.  The NOSA has been approved by the Open Source
 Initiative.  See the file COPYING at the top of the distribution
 directory tree for the complete NOSA document.

 THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
 KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
 LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
 SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
 A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
 THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
 DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.

 __END_LICENSE__
