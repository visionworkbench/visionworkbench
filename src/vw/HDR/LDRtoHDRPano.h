#ifndef __LDRTOHDRPANO_H__
#define __LDRTOHDRPANO_H__

// Vision Workbench
#include <vw/HDR/LDRtoHDR.h>
#include <vw/HDR/ExifView.h>

// STL
#include <list>

namespace vw { namespace HDR {

const int SECTION_SIZE = 256; // For final version, change to something like 2k

const int NUM_PAIRS_PANO = 5000;

struct Rect {
  int x;
  int y;
  int width;
  int height;

  Rect() : x(0), y(0), width(0), height(0) { }
  Rect(int x, int y, int w, int h) : x(x), y(y), width(w), height(h) { }
};

/* Essentially a placeholder for whatever ImageView(Ref) business is eventually used as inputs.
 * We need offsets and dimensions so we can figure out which images overlap without actually
 * loading them (memory inefficient for very large datasets).
 */
struct PanoImage {
  char * filename;
  Rect rect;
  double brightness_multiplier;

  PanoImage(char * name, Rect r, double b_m) : filename(name), rect(r), brightness_multiplier(b_m) { }
};

struct IntensityPair {
  PixelRGBA<double> I_1;
  PixelRGBA<double> I_2;
  double ratio;

  IntensityPair(PixelRGBA<double> i1, PixelRGBA<double> i2, double r) : I_1(i1), I_2(i2), ratio(r) { }
};

struct PanoLoadedImage {
  ImageView<PixelRGBA<double> > image;
  Rect rect;
  double Bv; // brightness value

  PanoLoadedImage(ImageView<PixelRGB<double> > i, Rect r, double b) : image(i), rect(r), Bv(b) { }

  bool contains(int x, int y) {
    return ((x >= rect.x) && (x < rect.x + rect.width) && (y >= rect.y) && (y < rect.y + rect.height));
  }

  bool is_defined(int x, int y) {
    printf("               Alpha: %f\n", image(x - rect.x, y - rect.y).a());
    return (contains(x, y) && (image(x - rect.x, y - rect.y).a() > 0));
  }

  PixelRGBA<double> get_pixel(int x, int y) {
    VW_ASSERT(contains(x, y), ArgumentErr() << "Point (" << x << "," << y << ") outside image bounds.");
    return image(x - rect.x, y - rect.y);
  }
};

void generate_intensity_pairs_pano(vector<PanoImage> &images, int pano_width, int pano_height, vector<double> const &brightness_values, vector<Matrix<double> > &pairs);

void process_ldr_images_pano(vector<PanoImage> &images, vector<double> const &brightness_values, int pano_width, int pano_height, vector<Vector<double> > &ret_curves) {
  vector<Matrix<double> > pairs(3);
  vector<Vector<double> > curves(3);

  // Do sampling
  generate_intensity_pairs_pano(images, pano_width, pano_height, brightness_values, pairs);

  // Compute camera response curves
  for (int i = 0; i < 3; i++) {
    estimate_camera_curve(pairs[i], curves[i], RESPONSE_POLYNOMIAL_ORDER);
  }

  ret_curves = curves;

  // Compute brightness_multipliers
  double min_bv = *(min_element(brightness_values.begin(), brightness_values.end()));
  for (int i = 0; i < images.size(); i++) {
    images[i].brightness_multiplier = pow(2.0, (brightness_values[i] - min_bv) * 0.5);
  }
}

void process_ldr_images_pano(vector<PanoImage> &images, int pano_width, int pano_height, vector<Vector<double> > &ret_curves) {
  vector<double> brightness_values(images.size());
  ExifView exif;

  // Get brightness values from Exif
  for (int i = 0; i < images.size(); i++) {
    VW_ASSERT(exif.load_exif(images[i].filename), IOErr() << "File " << images[i].filename << " contains no Exif data.");
    brightness_values[i] = exif.get_brightness_value();
    cout << "Brightness_value = " << brightness_values[i] << "\n";
  }

  process_ldr_images_pano(images, brightness_values, pano_width, pano_height, ret_curves);
}

void print_rect(Rect r) {
  printf("(%i, %i, %i, %i)\n", r.x, r.y, r.width, r.height);
}

Rect get_overlap(Rect r1, Rect r2) {
  Rect overlap;
  overlap.x = max(r1.x, r2.x);
  overlap.width = min(r1.x + r1.width, r2.x + r2.width) - overlap.x;
  overlap.y = max(r1.y, r2.y);
  overlap.height = min(r1.y + r1.height, r2.y + r2.height) - overlap.y;
  if ((overlap.width <= 0) || (overlap.height <= 0))
    return Rect(0, 0, 0, 0);
  return overlap;
}

bool is_overlap(Rect overlap) {
  return ((overlap.width > 0) && (overlap.height > 0));
}

bool is_overlap(Rect r1, Rect r2) {
  return is_overlap(get_overlap(r1, r2));
}

void generate_intensity_pairs_pano(vector<PanoImage> &images, int pano_width, int pano_height, vector<double> const &brightness_values, vector<Matrix<double> > &pairs) {
  typedef ImageView<PixelRGBA<double> > Image;

  int num_pairs_section = (int)((double)NUM_PAIRS_PANO / (ceil((double)pano_width / (double)SECTION_SIZE) * ceil((double)pano_height / (double)SECTION_SIZE)));
  printf("Approximate number of pairs per section: %i\n", num_pairs_section);
  list<IntensityPair> i_pairs;
  Rect section;

  srand(time(0)); // Initialize random number generator

  // Operate locally on SECTION_SIZE x SECTION_SIZE squares of the panorama
  for (section.x = 0; section.x < pano_width; section.x += SECTION_SIZE) {
    section.width = min(SECTION_SIZE, pano_width - section.x); // edge case
    for (section.y = 0; section.y < pano_height; section.y += SECTION_SIZE) {
      section.height = min(SECTION_SIZE, pano_height - section.y); // edge case
      
      // Load images that overlap with current section
      list<PanoLoadedImage> loaded;
      for (int i = 0; i < images.size(); i++) {
	if (is_overlap(section, images[i].rect)) {
	  Image image;
	  read_image(image, images[i].filename);
	  loaded.push_back(PanoLoadedImage(image, images[i].rect, brightness_values[i]));
	}
      }
      printf("      Loaded: %i\n", loaded.size());
      if (loaded.size() < 2) continue; // Need at least two images

      // Attempt to sample num_pairs_section points. We give up if it's taking too long.
      int attempts = 0;
      int samples = 0;
      int max_samples = min(num_pairs_section, section.width * section.height / 4);
      int max_attempts = max_samples * 8;
      while ((samples < max_samples) && (attempts < max_attempts)) {
	// Generate random indices
	int rand_x = section.x + GET_RAND(section.width);
	int rand_y = section.y + GET_RAND(section.height);
	attempts++;

	// Determine which images have pixel data at that point
	vector<PanoLoadedImage> defined;
	for (list<PanoLoadedImage>::iterator i = loaded.begin(); i != loaded.end(); i++) {
	  if ((*i).is_defined(rand_x, rand_y)) defined.push_back(*i);
	}

	int n_images = defined.size();
	printf("            Num images defined at (%i,%i): %i\n", rand_x, rand_y, n_images);
	if (n_images > 1) {
	  // Pick two distinct images to sample from
	  int id1 = GET_RAND(n_images);
	  int id2;
	  while (true) {
	    id2 = GET_RAND(n_images);
	    if (id1 != id2) break;
	  }

	  // Sample both images at those indices
	  PixelRGBA<double> I_1 = defined[id1].get_pixel(rand_x, rand_y);
	  PixelRGBA<double> I_2 = defined[id2].get_pixel(rand_x, rand_y);
	  // Calculate luminances (to prune out over- / under-exposed pixels)
	  PixelGray<double> gray_1(I_1);
	  PixelGray<double> gray_2(I_2);
	  double L_1 = gray_1.v();
	  double L_2 = gray_2.v();

	  // Add the sample pair to the list along with the exposure ratio
	  if ((L_1 > 0.01) && (L_1 < 0.99) && (L_2 > 0.1) && (L_2 < 0.9)) {
	    double ratio = pow(2.0, (defined[id1].Bv - defined[id2].Bv) * 0.5);
	    i_pairs.push_back(IntensityPair(I_1, I_2, ratio));
	    samples++;
	  }
	}
      }
    }
  }

  // Populate pairs matrices (for benefit of CameraCurve functions)
  // We could write an overloaded version of estimate_camera_curves instead.
  pairs.clear();
  for (int channel = 0; channel < 3; channel++) {
    Matrix<double> pair_mat(i_pairs.size(), 3);
    int i = 0;
    for (list<IntensityPair>::iterator iter = i_pairs.begin(); iter != i_pairs.end(); iter++) {
      pair_mat(i, 0) = (*iter).I_1[channel];
      pair_mat(i, 1) = (*iter).I_2[channel];
      pair_mat(i, 2) = (*iter).ratio;
      i++;
    }
    pairs.push_back(pair_mat);
  }
}

}} /* vw::HDR */

#endif  // __LDRtoHDRPANO_H__
