#!/usr/bin/env python
# __BEGIN_LICENSE__
# Copyright (C) 2006-2011 United States Government as represented by
# the Administrator of the National Aeronautics and Space Administration.
# All Rights Reserved.
# __END_LICENSE__


from __future__ import with_statement

import sys
import os
from subprocess import Popen, PIPE

class Graph(object):
    def __init__(self, seq):
        self.e = set()
        self.v1 = set()
        self.v2 = {}

        for line in seq:
            line = line.strip()
            if not line:
                continue
            left, right = line.split()
            self.e.add((left, right))
            self.v1.add(left)
            self.v1.add(right)

    def assign_depth(self):

        print >>sys.stderr, 'there are %i nodes' % len(self.v1)
        for v in self.leaves():
            self.v1.discard(v)
            self.v2[v] = 0

        last = 10000000000000000
        times = 0
        while self.v1:
            cur = len(self.v1)
            print >>sys.stderr, 'now there are %i nodes' % cur
            if cur == last:
                times += 1
                if times > 10:
                    raise Exception('ERROR: Cycle! Nodes: %s' % self.v1)
            else:
                last = cur
                times = 0

            kill = set()
            for v in self.v1:
                children = self.children_of(v)
                if set(children).issubset(self.v2.keys()):
                    self.v2[v] = max([self.v2[i] for i in children])+1
                    kill.add(v)
            self.v1.difference_update(kill)

        print >>sys.stderr, 'now there are 0 nodes'

    def leaves(self):
        r = set([i[1] for i in self.e])
        r.difference_update([i[0] for i in self.e])
        return r

    def children_of(self, v):
        return [i[1] for i in self.e if i[0] == v]

if __name__ == '__main__':

    args = len(sys.argv)

    infile  = None
    outfile = None

    if args == 3:
        outfile = sys.argv[1]
        infile  = sys.argv[2]
    elif args == 2:
        outfile = sys.argv[1]
    else:
        print >>sys.stderr, 'Usage: %s output-file' % sys.argv[0]
        sys.exit(-1)


    with file(outfile, 'w') as f:
        if infile is None:
            grep = Popen(['grep', '-rI', '#include.*<vw', os.path.abspath(os.path.dirname(sys.argv[0]) + '/../src/vw')], stdout=PIPE)
            sed  = Popen(['sed', '-e', 's#^.*src/\(vw[^:]\+\):.*<\([^>]\+\)>#\\1 \\2#'], stdin=grep.stdout, stdout=PIPE)
            filt = Popen(['grep', '-v', '\(.cc\|/tests\|/GPU\)'], stdin=sed.stdout, stdout=PIPE)
            outf = filt.communicate()[0]
        else:
            outf = file(infile, 'r').read()

        g = Graph(outf.split('\n'))
        g.assign_depth()

        print >>f, '\n'.join(map(lambda x: '%i\t%s' % (x[1], x[0]), sorted(g.v2.items(), key=lambda x: (x[1],x[0]))))
