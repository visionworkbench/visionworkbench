#!/usr/bin/perl
# __BEGIN_LICENSE__
# Copyright (C) 2006, 2007 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


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
    ".hh"   => "//",
    ".i"    => "//",
    ".m4"   => "dnl",
    ".mak"  => "#",
    ".pl"   => "#",
    ".py"   => "#",
    ".sh"   => "#",
    ".js"   => "//",
    ".tcc"  => "//",
);

# Read the license text from __DATA__ by default
my $f = \*DATA;
$f = $ARGV[0] if @ARGV > 0;

my @license = read_file($f);
my $shebang = '';

# process each line given on stdin
foreach my $filename (<>) {
    chomp $filename;

    # get the extension, and skip it if we don't know about it
    my (undef, undef, $ext) = fileparse($filename, qr/\.[^.]*/);

    unless (exists $comment{$ext}) {
        warn "Skipped $filename\n";
        next;
    }

    my $file = read_file($filename);

    $shebang = '';
    # Protect a shebang line
    if ($file =~ s/^(#!.*\n)//) {
        if (defined($1)) {
            $shebang = $1;
        }
    }

    # Remove a license header if it exists
    $file =~ s/^[^\n]*__BEGIN_LICENSE__.*?__END_LICENSE__[^\n]*$//ms;

    # Remove all blank files from the top of the file
    while ($file =~ s/^\s*\n//) {};

    # prepend the license text, prepending the comment string to each line.
    # Also, separate the license header from content by two blank lines
    $file = $shebang . $comment{$ext} . join($comment{$ext}, @license) . "\n\n" . $file;

    write_file($filename, $file);
}

__DATA__
 __BEGIN_LICENSE__
 Copyright (C) 2006, 2007 United States Government as represented by
 the Administrator of the National Aeronautics and Space Administration.
 All Rights Reserved.
 __END_LICENSE__
