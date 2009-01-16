#!/usr/bin/env python

import sys
import os

class Graph(object):
    def __init__(self, seq):
        self.e = set()
        self.v1 = set()
        self.v2 = {}

        for line in seq:
            left, right = line.split()
            self.e.add((left, right))
            self.v1.add(left)
            self.v1.add(right)

    def assign_depth(self):
        print >>sys.stderr, 'there are %i nodes' % len(self.v1)
        for v in self.leaves():
            self.v1.discard(v)
            self.v2[v] = 0

        while self.v1:
            print >>sys.stderr, 'now there are %i nodes' % len(self.v1)
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
    (inf, outf) = os.popen2(r"grep -r '#include.*<vw' " + os.path.abspath(os.path.dirname(sys.argv[0]) + '/../src/vw') + r" | sed -e 's#^.*src/\(vw[^:]\+\):.*<\([^>]\+\)>#\1 \2#' | grep -v '.cc ' | grep -v '/tests' | grep -v '/GPU'")
    inf.close()
    g = Graph(outf)
    outf.close()
    g.assign_depth()
    print '\n'.join(map(lambda x: '%i\t%s' % (x[1], x[0]), sorted(g.v2.items(), key=lambda x: (x[1],x[0]))))
