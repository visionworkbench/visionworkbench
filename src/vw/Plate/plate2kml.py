#!/usr/bin/env python2.7
import optparse
import sys, os
import math
sys.path.insert(1, '/Users/ted/local/src/libkml-1.2.0/build/lib/python2.7/site-packages') # tmo
import kmldom
import kmlengine
import subprocess
import json
import zipfile
import urlparse
import multiprocessing
import Queue # this is just imported for the Queue.Empty exception <shrug>


DEFAULT_MAX_KMZ_FEATURES = 256
DEFAULT_DRAW_ORDER_FACTOR = 10 # presently hardcoded: tiles will be rendererd at DRAW_ORDER_FACTOR * tile_level
DEFAULT_VW_BIN_PATH = '/Users/ted/local/bin/'
DEFAULT_REGIONATION_OFFSET = 5
DEFAULT_WORKERS = 6

def int_or_none(val):
    try:
        return int(val)
    except TypeError:
        return None

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
                maxx = minx + width
                maxy = miny + height

            if not allow_floats:
                for prop in initprops:
                    setattr(self, prop, int_or_none(locals()[prop]))
            else:
                for prop in initprops:
                    setattr(self, prop, locals()[prop])
    @property
    def width(self):
        return self.maxx - self.minx
    @width.setter
    def width(self, value):
        self.maxx = self.minx + value

    @property
    def height(self):
        return self.maxy - self.miny
    @height.setter
    def height(self, value):
        self.maxy = self.miny + value

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
        self.latlon_bbox = BBox(allow_floats=True) # initially empty, this decouples the region's degree-space bounds from it's tile-space bounds.

    # TODO: semantically, row and col should be minx and miny
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
        region = factory.CreateRegion()
        
        if self.level == 1:
          ## Top layerish. Layer 0 is never drawn because it's
          ## boundaries are illegal (Lat range too large). So this step
          ## is to make sure layer 1 can always be seen when zoomed out.
            #minlod, maxlod = (1,513)
            minlod, maxlod = (1,1024)
        # TODO: If we can figure out how to tell if this is the deepest plate level, it'd be great to set the maxlod to infinity (-1)
        #    elif lowest_overlay:
        #        # end of branch
        #        minlod, maxlod = (128 ,-1)
        else:
            #minlod, maxlod = (128, 531)
            #minlod, maxlod = (64, -1)
            #minlod, maxlod = (64, max((self.height, self.width)) * 2 + 1)
            minlod = int(math.sqrt(self.height * self.width)) * 32 - 4 #64
            maxlod = int(math.sqrt(self.height * self.width )) * 512 * 4 + 1
            #maxlod = -1
            if self.level == options.max_levels:
                maxlod = -1
        region.set_lod(create_lod(minlod, maxlod, factory))
        
        assert not self.latlon_bbox.is_null() # presumably, we already added some points to the latlon_bbox with latlon_bbox.expand()
        latlonalt_box = create_latlonalt_square(self.latlon_bbox.maxy, self.latlon_bbox.miny, self.latlon_bbox.minx, self.latlon_bbox.maxx, factory)
        region.set_latlonaltbox(latlonalt_box)
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

    def latlonbox(self):
        """
        lon_delta = 360.0 / (2**self.level)
        lat_delta = 180.0 / (2**self.level)
        west = -180 + self.col * lon_delta
        east = west + lon_delta
        north = 90 - lat_delta * self.row
        south = north - lat_delta
        """

        llbox = factory.CreateLatLonBox()
        llbox.set_west(self.west)
        llbox.set_south(self.south)
        llbox.set_east(self.east)
        llbox.set_north(self.north)
        return llbox

    def bbox(self):
        return BBox(
            self.col,
            self.row,
            1,
            1
        )    

    def kml_region(self, factory):
        region = factory.CreateRegion()
        
        if self.level == 1:
          ## Top layerish. Layer 0 is never drawn because it's
          ## boundaries are illegal (Lat range too large). So this step
          ## is to make sure layer 1 can always be seen when zoomed out.
            #minlod, maxlod = (1,513)
            minlod, maxlod = (1,1024)
        #    elif lowest_overlay:
        #        # end of branch
        #        minlod, maxlod = (128 ,-1)
        else:
            #minlod, maxlod = (128, 531)
            #minlod, maxlod = (64, -1)
            #minlod, maxlod = (64, 513)
            minlod, maxlod = (64, 513)
            if self.level == options.max_levels:
                maxlod = -1
        region.set_lod(create_lod(minlod, maxlod, factory))
        
        latlonalt_box = create_latlonalt_square(self.north, self.south, self.west, self.east, factory)
        region.set_latlonaltbox(latlonalt_box)
        return region


    def kml_overlay(self, factory):
        goverlay = factory.CreateGroundOverlay()
        #region = TileRegion(self.level, BBox(self.col, self.row, 1, 1)).kml_region()
        region = self.kml_region(factory)
        goverlay.set_region(region)
        goverlay.set_latlonbox(self.latlonbox())
        # TODO: maybe draw order shourld be a fcn of level
        goverlay.set_draworder(self.level * DEFAULT_DRAW_ORDER_FACTOR) 

        url = urlparse.urljoin( options.baseurl, '/'.join( str(t) for t in (self.level, self.col, self.row) ) + "."+self.filetype )
        icon = factory.CreateIcon()
        icon.set_href( str(url) )
        goverlay.set_icon( icon )

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
    while True:
        try:
            output = subprocess.check_output( args )
            break
        except subprocess.CalledProcessError, e:
            print "%s RETRYING" % str(e)
            continue
    response = json.loads(output)
    if response['ok']:
        tiles = [ Tile(**t) for t in response['result'] ]
    else:
        if response['error'] == "TOO DEEP":
            raise PlateDepthExceededException()
        else:
            raise PlateException(response['error'])
    print "FOUND %d" % len(tiles)
    return tiles

def create_lod(minpix, maxpix, factory):
    lod = factory.CreateLod()
    lod.set_minlodpixels(minpix)
    lod.set_maxlodpixels(maxpix)
    return lod

def create_latlonalt_square(north, south, west, east,factory):
    box = factory.CreateLatLonAltBox()
    box.set_south(south)
    box.set_north(north)
    box.set_west(west)
    box.set_east(east)
    return box
    
    
def overlay_for_tile(tile, factory):
    goverlay = factory.CreateGroundOverlay()
    region = TileRegion(tile.level, BBox(tile.col, tile.row, 1, 1)).kml_region()
    goverlay.set_region(region)
    goverlay.set_latlonbox(tile.latlonbox())
    goverlay.set_draworder(tile.level * DEFAULT_DRAW_ORDER_FACTOR) 

    url = urlparse.urljoin( options.baseurl, '/'.join( str(t) for t in (tile.level, tile.col, tile.row) ) + "."+tile.filetype )
    icon = factory.CreateIcon()
    icon.set_href( str(url) )
    goverlay.set_icon( icon )

    return goverlay

def make_netlink(url, region, factory, name=None):
    netlink = factory.CreateNetworkLink()
    if region:
        netlink.set_region(region.kml_region())
    link = factory.CreateLink()
    link.set_href(url)
    link.set_viewrefreshmode( kmldom.VIEWREFRESHMODE_ONREGION )
    netlink.set_link(link)
    if name:
        netlink.set_name(name)
    return netlink

def draw_level(levelno, plate, output_dir, prior_hit_regions, factory):
    level = PlateLevel(levelno)
    netlinks = []
    hit_regions = []
    
    region_queue = multiprocessing.JoinableQueue()
    output_queue = multiprocessing.Queue()

    print "Launching %d workers." % options.workers
    workers = []
    for i in range(options.workers):
        p = multiprocessing.Process( target=draw_region_worker, args=(region_queue, output_queue, plate, output_dir, factory) )
        workers.append(p)
        p.start()

    rcount = 0
    for region in level.regionate(restrict_to_regions=prior_hit_regions):
        try:
            #netlink = draw_region(region, plate, output_dir, factory)
            region_queue.put(region)
            rcount += 1
        except PlateDepthExceededException:
            print "Hit bottom!"
            break
    print "Searching %d regions at level %d" % (rcount, levelno)
    print "Waiting for join..."
    region_queue.join()
    print "Joined!"
    region_queue.close()

    for w in workers:
        w.terminate()

    while True:
        try:
            region = output_queue.get(False)
            netlink = make_netlink(os.path.basename(region.kmz_filename), region, factory, name=region.name())
            netlinks.append(netlink)
            hit_regions.append(region)
        except Queue.Empty:
            break
    output_queue.close()

    # clear the referenced list and replace contents with new hit_regions
    del prior_hit_regions[:]
    prior_hit_regions += hit_regions

    print "Tiles in %d regions at level %d" % (len(netlinks), levelno)

    if len(netlinks) == 0:
        return None
    else:
        outfilename = os.path.join(output_dir, "%02d"%level.level+".kml")
        with open(outfilename, 'w') as outfile:
            folder = factory.CreateFolder()
            folder.set_name("%02d" % level.level)
            for netlink in netlinks:
                folder.add_feature(netlink)
            kml = factory.CreateKml()
            kml.set_feature(folder)
            outfile.write( kmldom.SerializePretty( kml ) )

        #return make_netlink(outfilename, level.region, factory, name="%02d"%level.level)
        return make_netlink(os.path.basename(outfilename), None, factory, name="%02d"%level.level)


def draw_region(region, plate, output_dir, factory):
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
    for t in tiles:
        #goverlay = overlay_for_tile(t, factory)
        goverlay = t.kml_overlay(factory)
        goverlays.append(goverlay)

        region.latlon_bbox.expand(t.west,t.south)
        region.latlon_bbox.expand(t.east,t.north)

    folder = factory.CreateFolder()
    for goverlay in goverlays:
        folder.add_feature(goverlay)
    kml = factory.CreateKml()
    kml.set_feature(folder)
    kmlstr = kmldom.SerializePretty(kml)
    kmzfile = zipfile.ZipFile(os.path.join(output_dir, region.kmz_filename), 'w')
    kmzfile.writestr('doc.kml', kmlstr)
    kmzfile.close()

    return region

def draw_region_worker(region_queue, output_queue, plate, output_dir, factory):
    """
    Loop infinitely and try to draw regions.
    """
    while True:
        try:
            region = region_queue.get(True, 3)
        except Queue.Empty:
            #print "Region queue empty.  Retrying"
            continue
        try:
            result = draw_region(region, plate, output_dir, factory)
            if result:
                output_queue.put(result)
        finally:
            region_queue.task_done()

def draw_plate(platefile_url, output_dir):
    if not os.path.exists(output_dir):
        raise Exception("output dir does not exist: ", output_dir)
    global factory
    factory = kmldom.KmlFactory.GetFactory()
    level = 0
    keep_drilling = True
    netlinks = []
    hit_regions = []
    folder = factory.CreateFolder()
    kml = factory.CreateKml()
    if options.planet:
        kml.set_hint("target=%s" % options.planet)
    kml.set_feature(folder)
    while keep_drilling and level <= options.max_levels:
        print "Drawing level %d..." % level
        netlink = draw_level(level, platefile_url, output_dir, hit_regions, factory)
        keep_drilling = bool(netlink)
        if netlink: netlinks.append(netlink)

        for netlink in netlinks:
            folder.add_feature(netlink)
        print "Writing root.kml"
        with open(os.path.join(output_dir, 'root.kml'), 'w') as outfile:
            outfile.write(kmldom.SerializePretty(kml))

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


