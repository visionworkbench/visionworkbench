# __BEGIN_LICENSE__
#  Copyright (c) 2006-2012, United States Government as represented by the
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


# Core
import core
from core import ErrorMessage, WarningMessage, InfoMessage, DebugMessage, VerboseDebugMessage
from core import set_debug_level
from core import PythonProgressCallback as ProgressCallback, TerminalProgressCallback

# Math
import vwmath
from vwmath import BBox2i

# Pixel types
import pixel
from pixel import uint8, int16, uint16, float32
from pixel import PixelGray, PixelGrayA, PixelRGB, PixelRGBA, PixelHSV, PixelXYZ
from pixel import channel_type, pixel_type, pixel_format, channel_range

# Images
import image
from image import Image, isimage

# Pixel type casting for pixels and images
import pixelcast
from pixelcast import pixel_cast

# Edge Extension
# Note: The bounding-box form of edge_extend does not yet exist.
from image import ZeroEdgeExtension, ConstantEdgeExtension, PeriodicEdgeExtension, ReflectEdgeExtension
from image import edge_extend

# mage Math
import imagemath
from imagemath import sin, cos, tan

# Manipulation
import imagemanip
from imagemanip import transpose, rotate_180, rotate_90_cw, rotate_90_ccw, flip_vertical, flip_horizontal, crop, select_col, select_row, select_plane

# Algorithms
# Note: Grassfire is absent because it currently returns ImageView<int>, which we don't support.
import imagealgo
from imagealgo import fill, clamp, normalize, threshold, is_opaque

# Filters
# Note: We need better pythonic kernel representations.  Perhaps NumPy arrays?
import filter
from filter import convolution_filter, separable_convolution_filter, gaussian_filter, derivative_filter, laplacian_filter

# Transform
import transform
from transform import NearestPixelInterpolation, BilinearInterpolation, BicubicInterpolation
from transform import resample, resize, translate, rotate, homography

# File IO
import fileio
from fileio import DiskImageResource, DiskImageResourceJPEG, DiskImageResourcePNG, DiskImageView

# Quadtree
import qtree
from qtree import QuadTreeGenerator
import composite
from composite import ImageComposite

# Cartography
import cartography
from cartography import Datum, GeoReference, read_georeference, geotransform, PixelAsPoint, PixelAsArea

# Dispatch on the source resource type.
# Currently only filenames are supported.
def read_image( source, *args, **pargs ):
    return fileio.read_image( source, *args, **pargs )

# Dispatch on the destination resource type.
# Currently only filenames are supported.
def write_image( dest, *args, **pargs ):
    return fileio.write_image( dest, *args, **pargs )
