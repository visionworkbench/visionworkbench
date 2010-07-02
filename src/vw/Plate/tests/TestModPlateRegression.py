#!/usr/bin/env python

import unittest
import urllib2
import logging
import sys
import mimetypes

_magic_cookie = None
def sample_type(buf):
    import magic
    global _magic_cookie
    if _magic_cookie is None:
        _magic_cookie = magic.open(magic.MAGIC_MIME)
        _magic_cookie.load()
    typ = _magic_cookie.buffer(buf)
    if typ is None: return None
    return typ.split(';')[0]

mimetypes.add_type('application/x-wwt-automatic', '.auto', False)

def guess_type(filename):
    return mimetypes.guess_type(filename, strict=False)[0]

def parse_args():
    from optparse import OptionParser
    p = OptionParser()
    p.add_option('-q', action='store_const', help='Quiet', const=logging.WARNING, dest='loglevel', default=logging.INFO)
    p.add_option('-v', action='store_const', help='Loud',  const=logging.DEBUG, dest='loglevel')
    p.add_option('-w', dest='wtml', default='http://wwt.nasa.gov/static/wwt_mars.wtml', help='WTML to load [%default]')
    global opts, args
    (opts,args) = p.parse_args()

def setup_logging():
    #logging.basicConfig(level=opts.loglevel, format='%(levelname)-8s %(message)s', stream=sys.stdout)
    logging.basicConfig(level=opts.loglevel, format='%(levelname)-8s: %(message)s')
    global info, debug, warn
    info  = logging.info
    debug = logging.debug
    warn  = logging.warning


class GetReturn(object):
    def __init__(self, code, data, headers, time):
        setattr(self, 'code', code)
        setattr(self, 'data', data)
        setattr(self, 'headers', headers)
        setattr(self, 'time', time)

def get(url):
    return urllib2.urlopen(url)

def timed_get(url):
    debug('GET %s' % url)
    from time import time
    start = time()
    try:
        r = urllib2.urlopen(url)
    except urllib2.HTTPError, e:
        return GetReturn(code=e.code, data=None, headers=None, time=time()-start)

    return GetReturn(code=200, data=r.read(), headers=r.headers, time=time()-start)

def tile_get(url_template, level, col, row, dem=False):
    if dem:
        url = url_template.replace('{0}', str(level)).replace('{1}', str(col)).replace('{2}', str(row))
    else:
        url = url_template.replace('{1}', str(level)).replace('{2}', str(col)).replace('{3}', str(row))
    return timed_get(url)

'''
compare against cached proper root tile
Check for deep tile
    compare against cached proper tile

Check for unknown platefile
Check for nasty query strings
'''

class TestPlate(unittest.TestCase):
    def good_request(self, r):
        self.assertEqual(200, r.code)
        self.assertTrue(len(r.data) > 0)
        self.assertTrue(r.time < 1)

    def content_type(self, r, filename):
        g = guess_type(filename)
        if g != 'application/x-wwt-automatic':
            self.assertEqual(r.headers['Content-Type'], g)
        self.assertEqual(r.headers['Content-Type'], sample_type(r.data))

    def test_sanity(self):
        self.assertTrue(len(self.Url) > len(self.FileType))
        self.assertEqual(self.FileType, self.Url[-len(self.FileType):])
    def test_credits(self):
        self.assertTrue(len(self.Credits) > 0)
    def test_thumbnail(self):
        r = timed_get(self.ThumbnailUrl)
        self.good_request(r)
        self.content_type(r, self.ThumbnailUrl)
    def test_root(self):
        r = tile_get(self.Url, level=0, col=0, row=0)
        self.good_request(r)
        self.content_type(r, self.Url)
    def test_out_of_range(self):
        r = tile_get(self.Url, level=self.NumLevels, col=0, row=0)
        self.assertEqual(404, r.code)
        r = tile_get(self.Url, level=0, col=0, row=1)
        self.assertEqual(404, r.code)
        r = tile_get(self.Url, level=0, col=1, row=0)
        self.assertEqual(404, r.code)

if __name__ == '__main__':
    import new
    import xml.dom.minidom as DOM

    parse_args()
    setup_logging()

    suite = unittest.TestSuite()

    try:
        dom = DOM.parse(get(opts.wtml))
    except Exception, e:
        logging.critical('Failed to load and parse %s: %s' % (opts.wtml, e))
        sys.exit(-1)

    for imageset in dom.getElementsByTagName('ImageSet'):
        plate_info = dict([(k,imageset.attributes[k].value) for k in imageset.attributes.keys()])
        for elt in 'Credits', 'ThumbnailUrl':
            node = imageset.getElementsByTagName(elt)[0]
            if node.firstChild:
                plate_info[elt] = node.firstChild.wholeText

        # Empty elt != nonexistent elt, but this makes other code easier
        plate_info['__getattr__'] = lambda self, name: ''
        klass = new.classobj(str(plate_info['Name']), (TestPlate,), plate_info)
        tests = unittest.TestLoader().loadTestsFromTestCase(klass)
        suite.addTest(tests)

    sys.exit(0 if unittest.TextTestRunner(verbosity=1).run(suite).wasSuccessful() else 1)
