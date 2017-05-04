
'''
    Prepare a Sentinel-1 image and run a water detection tool on it.
    
    This tool takes care of two additional steps: 
      - Border junk removal.  This needs to be done prior to warping.
      - Image georegistration using gdalwarp.
    Water detection is run after these steps.
'''

import sys, os, optparse

# Append this file location to the system path to make sure that
#  the packaged GDAL executables are found.
os.environ['PATH'] += os.pathsep + os.path.dirname(os.path.realpath(__file__))

def copyGeoTiffInfo(inputGeoTiffPath, inputTiffPath, outputPath):
    '''Copies geo information from a GeoTiff to a plain Tiff'''

    # Extract geo information from the input geotiff
    geoInfoPath = inputGeoTiffPath + '.metadata'
    cmd = 'listgeo ' + inputGeoTiffPath +' > '+ geoInfoPath
    print cmd
    os.system(cmd)
    if not os.path.exists(geoInfoPath):
        raise Exception('Failed to extract metadata from input file!')
    
    # Paste the geo information into the output file
    cmd = 'geotifcp ' + inputTiffPath +' '+ outputPath +' -g '+ geoInfoPath
    print cmd
    os.system(cmd)
    
    # Clean up temporary file
    os.remove(geoInfoPath)

def main():
    
    # Get the input arguments from the command line
    usage  = "Usage: sentinel1_flood_detect.py <input path> <output path>\n"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--sensitivity", dest="sensitivity", default=None,
                      help="Decrease to detect more water, increase to detect less.")
    parser.add_option("--debug", dest="debug", default=False, action='store_true',
                      help="Enable debug output.")
    parser.add_option("--num-threads", dest="num_threads", default=None,
                      help="Number of threads to use.")
    parser.add_option("--tile-size", dest="tile_size", default=None,
                      help="Specify size of tile used for processing.")
    parser.add_option("--dem-path", dest="dem_path", default=None,
                      help="Process with this DEM file.")
    (options, args) = parser.parse_args()
    if len(args) < 2:
        print usage
        return
      
    # Parse input arguments
    input_path  = args[0]
    output_path = args[1]

    # Create output folder    
    output_folder = os.path.dirname(output_path)
    if output_folder == '':
        output_folder = '.'
    if not os.path.exists(output_folder):
      print 'Creating output folder: ' + output_folder
      os.mkdir(output_folder)

    # Set paths
    temp_path           = output_path + "_clean.tif"
    border_correct_path = output_path + "_border_correct.tif"
    ortho_path          = output_path + '_WARP.tif'

    # The ortho step is slow so skip it if the file exists
    if not os.path.exists(ortho_path):
    
        # Correct borders   
        cmd = 'clean_sentinel1_borders -o ' + temp_path +' '+ input_path
        print cmd
        os.system(cmd)
        if not os.path.exists(temp_path):
            raise Exception('Failed to clean image borders!')
        
        # Paste the metadata on to the border corrected image
        copyGeoTiffInfo(input_path, temp_path, border_correct_path)
        
        # Orthoproject the image   
        cmd = 'gdalwarp  -multi -r cubic ' + border_correct_path +' '+ ortho_path
        print cmd
        os.system(cmd)
        if not os.path.exists(ortho_path):
            raise Exception('Failed during gdalwarp step!')
    
   
    # Call the main processing function
    cmd = ('detect_water --mode sentinel1 -o ' + output_path +' '+ ortho_path)
    if options.dem_path:
      cmd += ' --dem-path ' + options.dem_path
    if options.debug:
      cmd += ' --debug'
    if options.sensitivity:
        cmd += ' --sensitivity ' + options.sensitivity
    if options.tile_size:
        cmd += ' --tile-size ' + options.tile_size
    if options.num_threads:
        cmd += ' --num-threads ' + options.num_threads
    print cmd
    os.system(cmd)
    if not os.path.exists(output_path):
        raise Exception('Failed during water detection step!')
    
    
    # Clean up intermediate files
    #os.remove(temp_path)
    #os.remove(border_correct_path)
    #os.remove(ortho_path)
    
    print 'Finished generating output file: ' + output_path

# When run from the command line, execute the main() function.
if __name__ == "__main__":
    sys.exit(main())
