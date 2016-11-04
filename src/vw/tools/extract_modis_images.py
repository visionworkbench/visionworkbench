#!/usr/bin/env python
# __BEGIN_LICENSE__
#  Copyright (c) 2006-2013, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NASA Vision Workbench is licensed under the Apache License,
#  Version 2.0 (the "License"); you may not use this file except in
#  compliance with the License. You may obtain a copy of the License at
#  http:#www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# __END_LICENSE__


import os, sys, optparse, subprocess


DESIRED_CHANNELS = ['sur_refl_b01_1',
                    'sur_refl_b02_1',
                    'sur_refl_b03_1',
                    'sur_refl_b04_1',
                    'sur_refl_b05_1',
                    'sur_refl_b06_1',
                    'sur_refl_b07_1',
                    'QC_500m_1',
                    'sur_refl_b01_1',
                    'sur_refl_b02_1',
                    'QC_250m_1'#,
                    # -- We need these, but GDAL won't open them!!!!!!!!!!!!!
                    #'state_1km_1',
                    #'gflags_1',
                    # -- We don't need these
                    #obscov_500m_1
                    #iobs_res_1
                    #q_scan_1
                    #num_observations
                    #num_observations_1km
                    #num_observations_500m
                    #obscov_1
                   ]

# TODO: Only write high res images from high res input!

def find_dataset_names(path):
    '''Returns the list of desired datasets contained in a file'''

    # Load info about the file
    cmd = ['gdalinfo', path]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    textOutput, err = p.communicate()

    #print textOutput
    datasets = []
    sizes    = []
    for line in textOutput.split():
        if not ('SUBDATASET' in line):
            continue
        if ('NAME' in line):
            start = line.find('=') + 1
            datasets.append(line[start:])
        else:
            start = line.find('[') + 1
            stop  = line.find('x', start)
            size  = int(line[start:stop])
            sizes.append(size)
        #print line[start:]

    
    if len(datasets) != len(sizes):
      raise Exception('Error: dataset sizes mismatch!')
    
    return zip(datasets, sizes)


def prune_datasets(datasets):
    '''Remove duplicate datasets and undesired datasets.'''
    
    outputs = []
    for d in datasets:
        name = d[0]
        size = d[1]

        # Check if the name is on the desired channel list
        found = False
        for c in DESIRED_CHANNELS:
            if c in name:
                found = True
                break
        if not found:
            continue
        
        # Look for duplicates
        keep = True
        for dO in datasets:
            nameO = dO[0]
            sizeO = dO[1]
            
            if ((name == nameO) and (sizeO > size)):
                keep = False
                break               
        if not keep:
            continue

        outputs.append(d)

    return outputs
    
    

def extract_datasets(path, datasets, options):
    '''Extract the desired datasets to the output folder'''
    
    for d in datasets:

        text = d[0]
        name = text[text.rfind(':')+1:]
        outputPath = options.prefix + name + '.tif'
        if (options.overwrite or (not os.path.exists(outputPath))):
            cmd = "gdal_translate -of GTiff '" + text + "' " + outputPath
            print cmd
            os.system(cmd)


def main():

    #try:
    usage = "usage: extract_modis_images.py [options]\n"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("--prefix", dest="prefix", default="",
                      help="Output prefix to use.")
    parser.add_option("--overwrite", dest="overwrite", default=False, action='store_true',
                      help="Overwrite existing output files.")

    (options, args) = parser.parse_args()

    if not args: 
        parser.error("need .input files")

    # Handle input arguments
    input_paths = args
    output_dir = os.path.dirname(options.prefix)
    if not os.path.exists(output_dir):
        os.system('mkdir -p ' + output_dir)

    print 'Starting processing'

    # Loop through the input files
    for path in input_paths:
    
        # Extract all the subdatasets
        datasets = find_dataset_names(path)
    
        datasets = prune_datasets(datasets)
        
        extract_datasets(path, datasets, options)
    
    print 'Finished!'


    #except:
    #    print usage
    #    return -1




if __name__ == "__main__":
    sys.exit(main())
