#!/usr/bin/env python2.7
import optparse
import sys, os
sys.path.insert(1, '/Users/ted/local/src/libkml-1.2.0/build/lib/python2.7/site-packages') # tmo
import kmldom
import kmlengine
import subprocess
import json
import zipfile
import urlparse


DEFAULT_MAX_KMZ_FEATURES = 256
DEFAULT_DRAW_ORDER = 50 # presently hardcoded
DEFAULT_VW_BIN_PATH = '/Users/ted/local/bin/'

#class dotdict(dict):
#    """Syntax hack"""
#    __getattr__ = dict.__getitem__
#    __setattr__= dict.__setitem__
#    __delattr__= dict.__delitem__

def delegate_property(delegate_object, propname, readonly=False):
    """Delegate attribute access to a member object."""
    def getter(self):
        """Access is delegated to self.%s""" % delegate_object
        return getattr(getattr(self, delegate_object), propname)
    if readonly:
        return property(getter)
    def setter(self, value):
        setattr(getattr(self, delegate_object), propname, value)
    return property(getter, setter)

class BBox(object):
    def __init__(self, minx, miny, width, height):
        self.minx = minx
        self.miny = miny
        self.width = width
        self.height = height

    @property
    def corners(self):
        return ( (self.minx, self.miny), (self.minx + self.width, self.miny, self.height) )

    def contains(self, bbox):
        b1 = self.corners
        b2 = bbox.corners
        return b1[0][0] <= b2[0][0] and\
        b1[0][1] <= b2[0][1] and\
        b1[1][0] >= b2[1][0] and\
        b1[1][1] >= b2[1][1]

    def quarter(self):
        """Break a region in to 4 pieces"""
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
        self.bbox = bbox
    col = delegate_property('bbox','minx',readonly=True)
    row = delegate_property('bbox','miny',readonly=True)
    width = delegate_property('bbox','width',readonly=True)
    height = delegate_property('bbox','height',readonly=True)

    def name(self):
        return "L%02dR%dC%d" % (self.level, self.row, self.col)
    def __hash__(self):
        return (self.level, self.bbox.minx, self.bbox.miny, self.bbox.width, self.bbox.height).__hash__()

    def subdivide(self):
        for bbox in self.bbox.quarter():
            yield self.TileRegion(self.level, bbox)

    def kml_region(self):
        deg_delta = 360.0 / (1 << self.level)
        lon_west = -180.0 + self.col*deg_delta
        lat_south = 180.0 - deg_delta*(self.row+1)
        region = factory.CreateRegion()
        
        if self.level == 1:
          ## Top layerish. Layer 0 is never drawn because it's
          ## boundaries are illegal (Lat range too large). So this step
          ## is to make sure layer 1 can always be seen when zoomed out.
            minlod, maxlod = (1,513)
        #    elif lowest_overlay:
        #        # end of branch
        #        minlod, maxlod = (128 ,-1)
        else:
            minlod, maxlod = (128, 531)
        region.set_lod(create_lod(minlod, maxlod, factory))
        
        latlonalt_box = create_latlonalt_square(lat_south, lon_west, deg_delta, factory)
        region.set_latlonaltbox(latlonalt_box)
        return region

class Tile(object):
    def __init__(self, row, col, level):
        self.row = row
        self.col = col
        self.level = level

    def latlonbox(self):
        deg_delta = 360.0 / (2**self.level)
        west = -180 + self.col * deg_delta
        south = 180 - deg_delta*(self.row+1)
        east = west + deg_delta
        north = south + deg_delta

        llbox = factory.CreateLatLonBox()
        llbox.set_west(west)
        llbox.set_south(south)
        llbox.set_east(east)
        llbox.set_north(north)
        return llbox


def search_tiles_by_region(platefile_url, region):
    args = os.path.join(options.vw_bin_path,'tiles4region') + ' %s %d %d %d %d -l %d' % (platefile_url, region.col, region.row, region.width, region.height, region.level)
    output = subprocess.check_output(args.split(' '))
    output = output.replace("'", '"')
    tiles = [ Tile(**t) for t in json.loads(output) ]
    return tiles

def create_lod(minpix, maxpix, factory):
    lod = factory.CreateLod()
    lod.set_minlodpixels(minpix)
    lod.set_maxlodpixels(maxpix)
    return lod

def create_latlonalt_square(south, west, degree_size, factory):
    box = factory.CreateLatLonAltBox()
    box.set_south(south)
    box.set_north(south + degree_size)
    box.set_west(west)
    box.set_east(west+degree_size)
    return box
    
    
def overlay_for_tile(tile, factory, extension="jpg"):
    goverlay = factory.CreateGroundOverlay()
    region = TileRegion(tile.level, BBox(tile.col, tile.row, 1, 1)).kml_region()
    goverlay.set_region(region)
    goverlay.set_latlonbox(tile.latlonbox())
    # TODO: maybe draw order shourld be a fcn of level
    goverlay.set_draworder(DEFAULT_DRAW_ORDER) 

    url = urlparse.urljoin( options.baseurl, '/'.join( str(t) for t in (tile.level, tile.row, tile.col) ) + "."+extension )
    icon = factory.CreateIcon()
    icon.set_href( url )
    goverlay.set_icon( icon )

    return goverlay

def netlink_for_kmz(url, region, factory):
    netlink = factory.CreateNetworkLink()
    netlink.set_region(region.kml_region())
    link = factory.CreateLink()
    link.set_href(url)
    link.set_viewrefreshmode( kmldom.VIEWREFRESHMODE_ONREGION )
    netlink.set_link(link)
    return netlink

def draw_kml_region(region, plate, output_dir, factory):
    netlinks = []
    goverlays = []
    tiles = search_tiles_by_region(plate, region)
    if len(tiles) == 0:
        return None
    if len(tiles) > options.max_features:
        subregions = region.subdivide()
        for subregion in subregions:
            netlink = draw_kml_region(subregion, plate, output_dir,  factory)
            if netlink:
                netlinks.append(netlink)
    else:
        for t in tiles:
            goverlay = overlay_for_tile(t, factory)
            goverlays.append(goverlay)
            #subregion = TileRegion(t.level+1, BBox(t.col*2, t.row*2, t.col*2+1, t.row*2+1))
            subregion = TileRegion(t.level+1, BBox(t.col*2, t.row*2, 2, 2 ))
            netlinks.append( draw_kml_region(subregion, plate, output_dir, factory) )

    folder = factory.CreateFolder()
    for goverlay in goverlays:
        folder.add_feature(goverlay)
    for netlink in netlinks:
        folder.add_feature(netlink)
    kml = factory.CreateKml()
    kml.set_feature(folder)
    kmlstr = kmldom.SerializePretty(kml)
    kmz_filename = region.name() + ".kmz"
    kmzfile = zipfile.ZipFile(os.path.join(output_dir, kmz_filename), 'w')
    kmzfile.writestr('doc.kml', kmlstr)
    kmzfile.close()

    netlink = netlink_for_kmz(kmz_filename, region, factory)
    return netlink


def start_kml(platefile_url, output_dir):
    if not os.path.exists(output_dir):
        raise Exception("output dir does not exist: ", output_dir)
    global factory
    factory = kmldom.KmlFactory.GetFactory()
    initial_region = TileRegion(0, BBox(0,0,1,1))
    draw_kml_region(initial_region, platefile_url, output_dir, factory)

def main():
    usage = "%s <platefile> <output_dir>" % sys.argv[0]
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('--max-features', dest='max_features', help="Maximum number of features to stuff into one kmz file.", default=DEFAULT_MAX_KMZ_FEATURES)
    parser.add_option('--vw-bin', dest="vw_bin_path", help="Location of VW bin files", default=DEFAULT_VW_BIN_PATH)
    parser.add_option('--baseurl', dest='baseurl', help="mod_plate base URL", default="http://example.tld/path/99999")
    global options
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.print_help()
        sys.exit(1)  
    (platefile_url, output_dir) = args
    start_kml(platefile_url, output_dir)


if __name__ == "__main__":
    main()


