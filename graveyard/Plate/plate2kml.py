#!/usr/bin/env python
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


"""
plate2kml.py

A tool for generating a KML image pyramid based on a remote plate.
Dependencies:
    - Vision Workbench, built with Plate support
    - pykml, available from pypi: http://pypi.python.org/pypi/pykml
    - lxml: http://lxml.de/
"""

import optparse
import sys, os, time
import math
import subprocess
import json
import zipfile
import urlparse
import multiprocessing
import threading
import Queue # this is just imported for the Queue.Empty exception <shrug>

from lxml import etree
from pykml.factory import KML_ElementMaker as KML
from pykml.factory import ATOM_ElementMaker as ATOM
from pykml.factory import GX_ElementMaker as GX
import pykml.parser

TILE_MIN_LOD, TILE_MAX_LOD = (64, 513)
HARD_CODED_TILE_SIZE = 256

DEFAULT_MAX_KMZ_FEATURES = 256
DEFAULT_DRAW_ORDER_FACTOR = 10 # presently hardcoded: tiles will be rendererd at DRAW_ORDER_FACTOR * tile_level
#DEFAULT_VW_BIN_PATH = '/Users/ted/local/bin/'
DEFAULT_VW_BIN_PATH = '/big/local/visionworkbench/bin'
DEFAULT_REGIONATION_OFFSET = 5
DEFAULT_WORKERS = 6

def int_or_none(val):
    try:
        return int(val)
    except TypeError:
        return None


# This copy of the OGC KML schema is distributed with pykml.  There's no guaruntee it won't go out of date.
SCHEMA = pykml.parser.Schema("ogckml22.xsd") 
def assert_valid(element_tree):
    #if not SCHEMA.validate(element_tree):
    #    raise AssertionError("ElementTree fails validation:  %s" % str(SCHEMA.schema.error_log))
    try: 
        val = SCHEMA.schema.assertValid(element_tree)
    except:
        print etree.tostring(element_tree, pretty_print=True)
        raise
    return val

class BBox(object):
    """
    >>> boxa = BBox(0,0,4,4)
    >>> boxa.minx, boxa.miny, boxa.maxx, boxa.maxy, boxa.width, boxa.height
    (0, 0, 3, 3, 4, 4)
    
    >>> boxb = BBox(1,1,1,1)
    >>> boxb.minx, boxb.miny, boxb.maxx, boxb.maxy, boxb.width, boxb.height
    (1, 1, 1, 1, 1, 1)

    >>> boxa.contains(boxb)
    True
    >>> boxb.contains(boxa)
    False
    >>> boxb.is_null()
    False
    >>> BBox().is_null()
    True
    >>> bb = BBox()
    >>> bb.expand(0,0)
    >>> bb.is_null()
    False

    >>> boxa.expand(16,16)
    >>> boxa.expand(-16,-16)
    >>> boxa.minx, boxa.miny, boxa.maxx, boxa.maxy, boxa.width, boxa.height
    (-16, -16, 16, 16, 33, 33)
    >>> [b.corners for b in boxa.quarter()]
    [((-16, -16), (-1, -1)), ((-16, 0), (-1, 15)), ((0, -16), (15, -1)), ((0, 0), (15, 15))]
    
    
    """
    class BBoxTypeError(Exception):
        """Raise if a method is called that doesn't make sense for this kind of BBox"""
        pass

    def __init__(self, minx=None, miny=None, width=None, height=None, maxx=None, maxy=None, discrete=True):
        initprops =  ('minx','miny','maxx','maxy')
        args = locals()
        self.discrete = discrete
        if all( args[prop] == None for prop in  ('minx','miny','maxx','maxy','width','height') ):
            # init props are null.  initialize an empty bounding box
            for prop in initprops:
                setattr(self, prop, None)
        else:
            if (minx == None) or (miny == None):
                raise ValueError("Required: minimum x and y values")
            if not bool(width and height) ^ bool(maxx and maxy):
                raise ValueError("You must provide either width & height or maximum x & y values (not both)")

            if width != None or height != None:
                assert width != None and height != None and minx != None and miny != None
                assert (not maxx and not maxy)
                if not self.discrete:
                    raise BBoxTypeError("Width and height only make sense for discrete BBoxes")
                maxx = minx + width - 1
                maxy = miny + height - 1

            if self.discrete:
                for prop in initprops:
                    setattr(self, prop, int_or_none(locals()[prop]))
            else:
                for prop in initprops:
                    setattr(self, prop, locals()[prop])
    @property
    def width(self):
        if  self.discrete:
            return self.maxx - self.minx + 1
        else:
            return abs(self.maxx - self.minx)

    @property
    def height(self):
        if self.discrete:
            return self.maxy - self.miny + 1
        else:
            return abs(self.maxy - self.miny)

    def is_null(self):
        if any( getattr(self,p) == None for p in ('minx','miny','maxx','maxy') ):
            assert all( getattr(self,p) == None for p in ('minx','miny','maxx','maxy') )
            return True
        else:
            return False

    @property
    def corners(self):
        return ( (self.minx, self.miny), (self.maxx, self.maxy ))

    def contains(self, bbox):
        """
        Check whether the given bounding box is inside this one and return a boolean.
        """
        b1 = self.corners
        b2 = bbox.corners
        return b1[0][0] <= b2[0][0] and\
        b1[0][1] <= b2[0][1] and\
        b1[1][0] >= b2[1][0] and\
        b1[1][1] >= b2[1][1]

    def expand(self, x,y):
        """
        Expand the bounding box if the coordinates lie outside its bounds.
        If the given coordinates lie outside the bounds of the BBox, 
        expand the BBox to cover them.
        """
        if (self.minx == None) or (x < self.minx):
            self.minx = x
        if (self.maxx == None) or (x > self.maxx):
            self.maxx = x
        if (self.miny == None) or (y < self.miny):
            self.miny = y
        if (self.maxy == None) or (y > self.maxy):
            self.maxy = y



    def quarter(self):
        """(generator) Break out into 4 subregions, if possible"""
        if self.width == 1 and self.height == 1:
            yield self
            raise StopIteration
        elif self.width == 1 or self.height == 1 or self.width != self.height:
            # This assumption simplifies matters greatly
            raise NotImplementedError
        else:
            for x in (self.minx, self.minx + self.width / 2):
                for y in (self.miny, self.miny + self.height / 2):
                    yield BBox(x, y, self.width/2, self.height/2)
        

class TileRegion(object):
    def __init__(self, level, bbox):
        self.level = level
        if callable(bbox):
            self._bbox = bbox()
        else:
            self._bbox = bbox
        self.latlon_bbox = BBox(discrete=False) # initially empty, this decouples the region's degree-space bounds from it's tile-space bounds.
        self.tile_degree_size = None

    @property
    def minx(self):
        return self._bbox.minx
    @property
    def miny(self):
        return self._bbox.miny
    @property
    def maxx(self):
        return self._bbox.minx
    @property
    def maxy(self):
        return self._bbox.miny
    @property
    def width(self):
        return self._bbox.width
    @property
    def height(self):
        return self._bbox.height

    def name(self):
        return "L%02dR%dC%dW%dH%d" % (self.level, self.miny, self.minx, self.width, self.height)

    @property
    def kmz_filename(self):
        return self.name()+".kmz"

    def __hash__(self):
        return (self.level, self.minx, self.bbox.miny, self.bbox.width, self.bbox.height).__hash__()

    def tiles(self):
        """
        Iterate over the tiles contained by this region.
        Some of these may not be in the plate file.
        If you want to iterate over tiles actually in the plate, use search_tiles_by_region.
        """
        c = self.minx
        r = self.miny
        for i in range(self.width):
            for j in range(self.height):
                yield Tile(r+j, c+i, self.level)

    def subdivide(self):
        """
        Quarter this region.
        """
        for bbox in self.bbox.quarter():
            subregion = self.TileRegion(self.level, bbox)
            yield subregion

    def project_to_level(self, level):
        """
        Return a TileRegion representing the extent of this TileRegion,
        projected onto a different level of the tile pyramid.
        """
        
        level_delta = level - self.level
        if level_delta == 0:
            return self
        scale_factor = 2 ** level_delta
        proj_bbox = BBox(
            self.minx * scale_factor,
            self.miny * scale_factor,
            self.width * scale_factor,
            self.height * scale_factor
        )
        if level_delta < 0:
            # Ensure that region bounds are still integers
            for prop in (proj_bbox.minx, proj_bbox.miny, proj_bbox.width, proj_bbox.height):
                assert prop % 1 == 0

        return TileRegion( level, proj_bbox )

    def kml_region(self):
        region = KML.Region()
        
        if self.level == 1:
          ## Top layerish. Layer 0 is never drawn because it's
          ## boundaries are illegal (Lat range too large). So this step
          ## is to make sure layer 1 can always be seen when zoomed out.
            #minlod, maxlod = (1,513)
            minlod, maxlod = (1,1024)
        else:
            degree_width = abs(self.latlon_bbox.maxx - self.latlon_bbox.minx)
            degree_height = abs(self.latlon_bbox.maxy - self.latlon_bbox.miny)
            region_degree_size = math.sqrt(degree_width * degree_height)
            minlod = int(TILE_MIN_LOD / self.tile_degree_size * region_degree_size)
            maxlod = int(math.ceil(TILE_MAX_LOD / self.tile_degree_size * region_degree_size))


            # It isn't awesome that this relies on the command-line options to know when it's at the deepest level.
            # Perhaps there is a way to get this information from the Plate directly?
            if self.level == options.max_levels:
                maxlod = -1

        
        assert not self.latlon_bbox.is_null() # presumably, we already added some points to the latlon_bbox with latlon_bbox.expand()
        latlonalt_box = create_latlonalt_square(self.latlon_bbox.maxy, self.latlon_bbox.miny, self.latlon_bbox.minx, self.latlon_bbox.maxx)
        region.LatLonAltBox = latlonalt_box
        region.Lod = create_lod(minlod, maxlod)
        assert_valid(region)
        return region

class Tile(object):
    def __init__(self, row, col, level, west=None, east=None, north=None, south=None, filetype="png"):
        self.row = row
        self.col = col
        self.level = level

        self.west = west
        self.east = east
        self.north = north
        self.south = south
        self.filetype = filetype

        self.pixel_width = HARD_CODED_TILE_SIZE
        self.pixel_height = HARD_CODED_TILE_SIZE

    @property
    def degree_width(self):
        return abs(self.east - self.west)
    @property
    def degree_height(self):
        return abs(self.north - self.south)

    def latlonbox(self):
        """
        lon_delta = 360.0 / (2**self.level)
        lat_delta = 180.0 / (2**self.level)
        west = -180 + self.col * lon_delta
        east = west + lon_delta
        north = 90 - lat_delta * self.row
        south = north - lat_delta
        """

        llbox = KML.LatLonBox(
            KML.north(self.north),
            KML.south(self.south),
            KML.east(self.east),
            KML.west(self.west),
        )
        assert_valid(llbox)
        return llbox

    def bbox(self):
        return BBox(
            self.col,
            self.row,
            1,
            1
        )    

    def kml_region(self):
        
        if self.level == 1:
          ## Top layerish. Layer 0 is never drawn because it's
          ## boundaries are illegal (Lat range too large). So this step
          ## is to make sure layer 1 can always be seen when zoomed out.
            #minlod, maxlod = (1,513)
            minlod, maxlod = (1,1024)
        else:
            minlod, maxlod = (TILE_MIN_LOD, TILE_MAX_LOD)

            # Again, not awesome.  See comment for Tile's kml_region()
            if self.level == options.max_levels:
                maxlod = -1

        region = KML.Region()
        region.LatLonAltBox = create_latlonalt_square(self.north, self.south, self.west, self.east)
        region.Lod = create_lod(minlod, maxlod)
        
        assert_valid(region)
        return region


    def kml_overlay(self):
        url = urlparse.urljoin( options.baseurl, '/'.join( str(t) for t in (self.level, self.col, self.row) ) + "."+self.filetype )

        goverlay = KML.GroundOverlay()
        goverlay.Region = self.kml_region()
        goverlay.drawOrder = KML.drawOrder(self.level * DEFAULT_DRAW_ORDER_FACTOR)
        goverlay.Icon = KML.Icon(KML.href(url))
        goverlay.LatLonBox = self.latlonbox()


        assert_valid(goverlay)
        return goverlay
        
class PlateLevel(object): 
    def __init__(self, level):
        self.level = level
        self.region = TileRegion(level, BBox(0,0,2**level,2**level))

    def regionate(self, offset=-9999, restrict_to_regions=[]):
        """
        Generate regions for this plate level, based on some heuristics.
        If a list of restrict_to_regions is given, only produce regions contained by those regions.
        """
        if offset == -9999:
            offset = options.regionation_offset
        region_level = self.level - offset

        if region_level <= 0:
            yield self.region
        else:

            if restrict_to_regions:
                for reg in restrict_to_regions:
                    reg = reg.project_to_level(region_level)
                    for tile in reg.tiles():
                        yield TileRegion(region_level, tile.bbox).project_to_level(self.level)
            else:
                offset_region = TileRegion(region_level, BBox(
                    0,
                    0,
                    2**region_level,
                    2**region_level
                ))
                for tile in offset_region.tiles():
                    yield TileRegion(region_level, tile.bbox).project_to_level(self.level)

class PlateException(Exception):
    pass
class PlateDepthExceededException(PlateException):
    pass

def search_tiles_by_region(platefile_url, region):
    args = '%d %d %d %d -l %d' % (region.minx, region.miny, region.width, region.height, region.level)
    print "Searching for tiles at %s: " % region.name(),
    args = [os.path.join(options.vw_bin_path,'tiles4region'), platefile_url] + args.split(' ')
    #print " ".join(args)
    while True:
        try:
            #output = subprocess.check_output( args )
            p = subprocess.Popen( args , stdout=subprocess.PIPE)
            output, stdin = p.communicate()
            break
        except subprocess.CalledProcessError, e:
            print "%s RETRYING" % str(e)
            continue
    response = json.loads(output)
    if response['ok']:
        #tiles = [ Tile(**t) for t in response['result'] ]
        tiles = []
        for t in response['result']:
            t = dict((str(k), v) for k,v in t.items())
            tiles.append(Tile(**t))
    else:
        if response['error'] == "TOO DEEP":
            raise PlateDepthExceededException()
        else:
            raise PlateException(response['error'])
    print "FOUND %d" % len(tiles)
    return tiles

def create_lod(minpix, maxpix):
    lod = KML.Lod(
        KML.minLodPixels(minpix),
        KML.maxLodPixels(maxpix),
    )
    assert_valid(lod)
    return lod

def create_latlonalt_square(north, south, west, east):
    box = KML.LatLonAltBox(
        KML.north(north),
        KML.south(south),
        KML.east(east),
        KML.west(west),
    )
    #assert_valid(box)
    return box
    
    
def overlay_for_tile(tile):
    goverlay = KML.GroundOverlay()
    goverlay.Region = TileRegion(tile.level, BBox(tile.col, tile.row, 1, 1)).kml_region()
    goverlay.LatLonBox = tile.latlonbox()
    goverlay.drawOrder = KML.drawOrder(tile.level * DEFAULT_DRAW_ORDER_FACTOR)

    url = urlparse.urljoin( options.baseurl, '/'.join( str(t) for t in (tile.level, tile.col, tile.row) ) + "."+tile.filetype )
    icon = factory.CreateIcon()
    icon.set_href( str(url) )
    goverlay.set_icon( icon )

    assert_valid(goverlay)
    return goverlay

def make_netlink(url, region, name=None):
    netlink = KML.NetworkLink()
    if name:
        netlink.name = KML.name(name)
    if region:
        netlink.Region = region.kml_region()
    netlink.Link = KML.Link(
            KML.href(url),
            KML.viewRefreshMode('onRegion'),
        )

    assert_valid(netlink)
    return netlink

def draw_level(levelno, plate, output_dir, prior_hit_regions):
    level = PlateLevel(levelno)
    netlinks = []
    hit_regions = []

    region_queue = Queue.Queue()
    output_queue = Queue.Queue()
    die_event = threading.Event()

    print "Launching %d workers." % options.workers
    workers = []
    for i in range(options.workers):
        p = threading.Thread(target=draw_region_worker, args=(die_event, region_queue, output_queue, plate, output_dir))
        workers.append(p)
        p.start()

    rcount = 0
    for region in level.regionate(restrict_to_regions=prior_hit_regions):
        try:
            region_queue.put(region)
            rcount += 1
        except PlateDepthExceededException:
            print "Hit bottom!"
            break
    print "Searching %d regions at level %d" % (rcount, levelno)
    print "Waiting for join..."
    region_queue.join()
    print "Joined!"

    die_event.set()


    while True:
        try:
            time.sleep(0.1)
            region = output_queue.get(False)
            netlink = make_netlink(os.path.basename(region.kmz_filename), region, name=region.name())
            netlinks.append(netlink)
            hit_regions.append(region)
            output_queue.task_done()
        except Queue.Empty:
            break

    # clear the referenced list and replace contents with new hit_regions
    del prior_hit_regions[:]
    prior_hit_regions += hit_regions

    print "Tiles in %d regions at level %d" % (len(netlinks), levelno)

    if len(netlinks) == 0:
        return None
    else:
        outfilename = os.path.join(output_dir, "%02d"%level.level+".kml")
        with open(outfilename, 'w') as outfile:
            folder = KML.Folder(KML.name( "%02d" % level.level) )
            for netlink in netlinks:
                folder.append(netlink)
            kml = KML.kml(folder)
            assert_valid(kml)
            outfile.write( etree.tostring( kml, pretty_print=True ) )

        return make_netlink(os.path.basename(outfilename), None, name="%02d"%level.level)


def draw_region(region, plate, output_dir):
    """
    Pull a TileRegion from the region_queue.
    Render the TileRegion to KML and write it to disk.
    Create a libkml NetworkLink object for the file and put it on the netlink_queue.
    """

    goverlays = []

    tiles = search_tiles_by_region(plate, region)
    if len(tiles) == 0:
        return None
    #if len(tiles) > options.max_features:  # TODO: Implement (or remove) max features / netlinks

    t = tiles[0]
    region.tile_degree_size = math.sqrt(t.degree_width * t.degree_height) # this is actually consistent for the whole level.
    for t in tiles:
        goverlay = t.kml_overlay()
        goverlays.append(goverlay)

        region.latlon_bbox.expand(t.west,t.south)
        region.latlon_bbox.expand(t.east,t.north)

    folder = KML.Folder()
    for goverlay in goverlays:
        folder.append(goverlay)

    kml = KML.kml(folder)
    assert_valid(kml)

    kmlstr = etree.tostring(kml, pretty_print=True)
    kmzfile = zipfile.ZipFile(os.path.join(output_dir, region.kmz_filename), 'w')
    kmzfile.writestr('doc.kml', kmlstr)
    kmzfile.close()

    return region

def draw_region_worker(die_event, region_queue, output_queue, plate, output_dir):
    """
    Loop infinitely and try to draw regions.
    """
    while True:
        if die_event.is_set():
            break
        try:
            region = region_queue.get(True, 0.1)
        except Queue.Empty:
            #print "Region queue empty.  Retrying"
            continue
        try:
            result = draw_region(region, plate, output_dir)
            if result:
                output_queue.put(result)
        finally:
            time.sleep(0.01)
            region_queue.task_done()

def draw_plate(platefile_url, output_dir):
    if os.path.exists(output_dir):
        print "Output directory already exists: %s" % output_dir
        print "Replace? (y/n)",
        answer = raw_input()
        if answer.lower() != 'y':
            exit(0)
    else:
        print "Creating %s" % output_dir
        os.makedirs(output_dir)

    level = 0
    keep_drilling = True
    netlinks = []
    hit_regions = []
    folder = KML.Folder()
    kml = KML.kml(folder)
    if options.planet:
        kml.set('hint', "target=%s" % options.planet)
    while keep_drilling and level <= options.max_levels:
        print "Drawing level %d..." % level
        netlink = draw_level(level, platefile_url, output_dir, hit_regions)
        keep_drilling = bool(netlink)
        if netlink: netlinks.append(netlink)

        for netlink in netlinks:
            folder.append(netlink)
        print "Writing root.kml"
        with open(os.path.join(output_dir, 'root.kml'), 'w') as outfile:
            outfile.write(etree.tostring( kml, pretty_print=True ) )

        level += 1

    print "Done."
    sys.exit(0)

        
global options
def main():
    usage = "%s <platefile> <output_dir>" % sys.argv[0]
    parser = optparse.OptionParser(usage=usage)
    # TODO: Implement max_features
    #parser.add_option('--max-features', dest='max_features', help="Maximum number of features to stuff into one kmz file.", default=DEFAULT_MAX_KMZ_FEATURES)
    parser.add_option('--vw-bin', dest="vw_bin_path", help="Location of VW bin files", default=DEFAULT_VW_BIN_PATH)
    parser.add_option('--baseurl', dest='baseurl', help="mod_plate base URL", default="http://example.tld/path/99999")
    parser.add_option('-r','--regionation-offset', dest='regionation_offset', default=DEFAULT_REGIONATION_OFFSET)
    parser.add_option('-l','--max-level', dest='max_levels', type='int', default=9999, help="Stop drawing after this level. **It is desirable to speficy this, even if you're going to max resolution, because we use it to set the LOD for the last level of the pyramid.**")
    parser.add_option('--planet', dest='planet', help="Tell GE to load this in a particular planet mode (Earth, Moon, Mars)", default="mars")
    parser.add_option('--workers', '-w', dest='workers', type='int', help="Number of worker subprocesses to spawn", default=DEFAULT_WORKERS)
    #parser.add_option('--resume', dest='resume', action="store_true", default=False, help="Don't overwrite KMZ files if they already exist.")
    global options
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit(1)  
    (platefile_url, output_dir) = args
    draw_plate(platefile_url, output_dir)


if __name__ == "__main__":
    main()


