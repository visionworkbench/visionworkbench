
import sys
import os


'''
    Prepare a Sentinel-1 image and run a water detection tool on it.
'''


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
    if len(sys.argv) < 3:
        print 'Usage: sentinel1_flood_detect.py <input path> <output path>'
        return
    
    # We need the first two arguments but all others will be passed directly to the main c++ tool.
    input_path  = sys.argv[1]
    output_path = sys.argv[2]

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
        cmd = '/home/smcmich1/repo/visionworkbench/src/vw/tools/clean_sentinel1_borders -o ' + temp_path +' '+ input_path
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
    
    # TODO: Handle thread and tile size inputs?
    
    # Call the main processing function
    cmd = '/home/smcmich1/repo/visionworkbench/src/vw/tools/detect_water --mode sentinel1 -o ' + output_path +' '+ ortho_path
    print cmd
    os.system(cmd)
    if not os.path.exists(output_path):
        raise Exception('Failed during water detection step!')
    
    
    # Clean up intermediate files
    os.remove(temp_path)
    os.remove(border_correct_path)
    os.remove(ortho_path)
    
    
    print 'Finished generating output file: ' + output_path

if __name__ == "__main__":
    sys.exit(main())
