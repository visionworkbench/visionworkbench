#!/usr/bin/env python

import os, sys

if not os.path.exists("image2qtree"):
    os.mkdir("image2qtree")
os.chdir("image2qtree")
os.system("../geotif-generate")

# Create the image2qtree results
files = os.listdir(".")
for i in files:
    os.system("image2qtree -m kml %s" % i);

# Create master kml
mkml = open('master.kml','w')

mkml.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
mkml.write("<kml xmlns=\"http://earth.google.com/kml/2.2\">\n")
mkml.write("<Document>\n")
mkml.write("  <name>ListStyle radiofolder</name>\n")
mkml.write("  <Folder>\n")
mkml.write("    <name>Image2qtree results</name>\n")
mkml.write("    <Style>\n")
mkml.write("      <ListStyle>\n")
mkml.write("        <listItemType>radioFolder</listItemType>\n")
mkml.write("      </ListStyle>\n")
mkml.write("    </Style>\n")

for i in files:
    prefix = i[:-4]
    mkml.write("    <NetworkLink>\n")
    mkml.write("      <name>%s</name><refreshVisibility>1</refreshVisibility>\n" % i)
    mkml.write("      <Link><href>%s/%s.kml</href></Link>\n" % (prefix, prefix))
    mkml.write("    </NetworkLink>\n")

mkml.write("  </Folder>\n")
mkml.write("</Document>\n")
mkml.write("</kml>\n")
