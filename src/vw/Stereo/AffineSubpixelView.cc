#include <vw/Image/Interpolation.h>
#include <vw/Stereo/AffineSubpixelView.h>

using namespace vw;
using namespace vw::stereo;

inline double huber_robust_coefficient (double delta_norm, double b) {
  if (delta_norm < b)
    return delta_norm*delta_norm;
  else
    return 2*b*delta_norm - b*b;
}

inline double cauchy_robust_coefficient (double delta_norm, double b) {
  double b_sqr = b*b;
  return b_sqr*log(1+delta_norm*delta_norm/b_sqr);
}

inline double blake_zisserman_robust_coefficient (double delta_norm, double b) {
  return -log(exp(-(delta_norm*delta_norm) ) + b);
}


ImageView<float> AffineSubpixelView::compute_gaussian_weight_image(int kern_width, int kern_height) const {
  
  int center_pix_x = kern_width/2;
  int center_pix_y = kern_height/2;
  int two_sigma_sqr = 2*pow(kern_width/3,2);
  
  ImageView<float> weight(kern_width, kern_height);
  for (int j = 0; j < kern_height; ++j) {
    for (int i = 0; i < kern_width; ++i ) {
      weight(i,j) = exp(-1 * (pow(i-center_pix_x,2) + pow(j-center_pix_y,2)) / two_sigma_sqr);
    }
  }
  return weight;
}



AffineSubpixelView::pixel_type AffineSubpixelView::get_cached_pixel(double x, double y, int32 p) const {
  
  // This is the maximum number of pixels that the solution can be
  // adjusted by affine subpixel refinement.
  static const int AFFINE_SUBPIXEL_MAX_TRANSLATION = 5;

  // Skip pixels along the border
  if (x < m_left_bbox.min().x() || x >= m_left_bbox.max().x() || y < m_left_bbox.min().y() || y >= m_left_bbox.max().y()) 
    return PixelDisparity<float>();
  
  // Skip over pixels for which we have no initial disparity estimate
  if (m_disparity_map(x,y).missing())
    return PixelDisparity<float>();

  BBox2i current_window(x-m_kern_width/2, y-m_kern_height/2, m_kern_width, m_kern_height);
  Vector2 base_offset( -m_disparity_map(x,y).h() , -m_disparity_map(x,y).v() );          
        
  // Initialize our affine transform with the identity.  The
  // entries of d are laid out in row major order:
  // 
  //   | d(0) d(1) d(2) | 
  //   | d(3) d(4) d(5) |
  //   |  0    0    1   |
  //
  Vector<float,6> d;
  d(0) = 1.0;
  d(4) = 1.0;
  Vector2 offset;

  // Compute the derivative image patches
  CropView<CropView<ImageView<float> > > left_image_patch = crop(*m_left_cached_log_image, current_window);
  CropView<CropView<ImageView<float> > > I_x = crop(*m_x_deriv, current_window);
  CropView<CropView<ImageView<float> > > I_y = crop(*m_y_deriv, current_window);
        
  // Compute the base weight image
  int good_pixels = adjust_weight_image(m_weight, crop(m_disparity_map, current_window), m_weight_template);
        
  // Skip over pixels for which there are very few good matches
  // in the neighborhood.
  if (good_pixels < m_weight_threshold) 
    return PixelDisparity<float>();
                
  // Iterate until a solution is found or the max number of
  // iterations is reached.
  for (unsigned iter = 0; iter < 10; ++iter) {
    offset(0) = d[2];
    offset(1) = d[5];
        
    // First we check to see if our current subpixel translation
    // is less than one half of the window width.  If not, then
    // we are probably having trouble converging and we abort
    // this pixel!!
    if (norm_2(offset) > AFFINE_SUBPIXEL_MAX_TRANSLATION) 
      break;
          
    InterpolationView<EdgeExtensionView<CropView<ImageView<float> >, ZeroEdgeExtension>, BilinearInterpolation> right_interp_image =
      interpolate(*m_right_cached_log_image, BilinearInterpolation(), ZeroEdgeExtension());
        
    float x_base = x + m_disparity_map(x,y).h();
    float y_base = y + m_disparity_map(x,y).v();
    //          float error_total = 0;

    Matrix<float,6,6> rhs;
    Vector<float,6> lhs;
    for (int jj = -m_kern_height/2; jj <= m_kern_height/2; ++jj) {
      for (int ii = -m_kern_width/2; ii <= m_kern_width/2; ++ii) {
        int i = ii + m_kern_width/2;
        int j = jj + m_kern_height/2;

        // First we compute the pixel offset for the right image
        // and the error for the current pixel.
        float xx = x_base + d[0] * ii + d[1] * jj + offset(0);
        float yy = y_base + d[3] * ii + d[4] * jj + offset(1);
        float I_e_val = right_interp_image(xx,yy) - left_image_patch(i,j) + 1e-16; 
        //              error_total += pow(I_e_val,2);

        // Apply the robust cost function.  We use a huber
        // function to gently remove outliers for small errors,
        // but we set a hard limit a 5 times the cost threshold
        // to remove major (salt&pepper) noise.
        double thresh = 1e-3;

        // Cauchy seems to work well with thresh ~= 1e-4
        float robust_weight = sqrt(cauchy_robust_coefficient(fabs(I_e_val),thresh))/fabs(I_e_val);
        
        // Huber seems to work well with thresh >= 1e-5
        //        float robust_weight = sqrt(huber_robust_coefficient(fabs(I_e_val),thresh))/fabs(I_e_val);

        // Disable robust cost function altogether
        //        float robust_weight = 1;

        // We combine the error value with the derivative and
        // add this to the update equation.
        float I_x_val = robust_weight * m_weight(i,j) * I_x(i,j);
        float I_y_val = robust_weight * m_weight(i,j) * I_y(i,j);
        float I_x_sqr = robust_weight * m_weight(i,j) * I_x(i,j) * I_x(i,j);
        float I_y_sqr = robust_weight * m_weight(i,j) * I_y(i,j) * I_y(i,j);
        float I_x_I_y = robust_weight * m_weight(i,j) * I_x(i,j) * I_y(i,j);

        // Left hand side
        lhs(0) += ii * I_x_val * I_e_val;
        lhs(1) += jj * I_x_val * I_e_val;
        lhs(2) +=      I_x_val * I_e_val;
        lhs(3) += ii * I_y_val * I_e_val;
        lhs(4) += jj * I_y_val * I_e_val;
        lhs(5) +=      I_y_val * I_e_val;
              
        // Right Hand Side UL
        rhs(0,0) += ii*ii * I_x_sqr;
        rhs(0,1) += ii*jj * I_x_sqr;
        rhs(0,2) += ii    * I_x_sqr;
        rhs(1,1) += jj*jj * I_x_sqr;
        rhs(1,2) += jj    * I_x_sqr;
        rhs(2,2) +=         I_x_sqr;
            
        // Right Hand Side UR
        rhs(0,3) += ii*ii * I_x_I_y;
        rhs(0,4) += ii*jj * I_x_I_y;
        rhs(0,5) += ii    * I_x_I_y;
        rhs(1,4) += jj*jj * I_x_I_y;
        rhs(1,5) += jj    * I_x_I_y;
        rhs(2,5) +=         I_x_I_y;
            
        // Right Hand Side LR
        rhs(3,3) += ii*ii * I_y_sqr;
        rhs(3,4) += ii*jj * I_y_sqr;
        rhs(3,5) += ii    * I_y_sqr;
        rhs(4,4) += jj*jj * I_y_sqr;
        rhs(4,5) += jj    * I_y_sqr;
        rhs(5,5) +=         I_y_sqr;
      }
    }

    lhs *= -1;

    // Fill in symmetric entries
    rhs(1,0) = rhs(0,1);
    rhs(2,0) = rhs(0,2);
    rhs(2,1) = rhs(1,2);
    rhs(1,3) = rhs(0,4);
    rhs(2,3) = rhs(0,5);
    rhs(2,4) = rhs(1,5);
    rhs(3,0) = rhs(0,3);
    rhs(3,1) = rhs(1,3);
    rhs(3,2) = rhs(2,3);
    rhs(4,0) = rhs(0,4);
    rhs(4,1) = rhs(1,4);
    rhs(4,2) = rhs(2,4);
    rhs(4,3) = rhs(3,4);
    rhs(5,0) = rhs(0,5);
    rhs(5,1) = rhs(1,5);
    rhs(5,2) = rhs(2,5);
    rhs(5,3) = rhs(3,5);
    rhs(5,4) = rhs(4,5);

    //           {          
    //             ImageView<float> right_image_patch(kern_width, kern_height);
    //             for (int jj = -kern_height/2; jj <= kern_height/2; ++jj) {
    //               for (int ii = -kern_width/2; ii <= kern_width/2; ++ii) {
    //                 float xx = x_base + d[0] * ii + d[1] * jj + offset(0);
    //                 float yy = y_base + d[3] * ii + d[4] * jj + offset(1);
    //                 right_image_patch(ii+kern_width/2, jj+kern_width/2) = right_interp_image(xx,yy);
    //               }
    //             }
    //             std::ostringstream ostr;
    //             ostr << x << "_" << y << "-" << iter;
    //             write_image("small/left-"+ostr.str()+".tif", left_image_patch);
    //             write_image("small/right-"+ostr.str()+".tif", right_image_patch);
    //             write_image("small/weight-"+ostr.str()+".tif", w);
    //           }
        

    // Solves lhs = rhs * x, and stores the result in-place in lhs.
    Matrix<double,6,6> pre_rhs = rhs;
    Vector<double,6> pre_lhs = lhs;
    try { 
      solve_symmetric_nocopy(rhs,lhs);
    } catch (ArgumentErr &e) {
      std::cout << "Error @ " << x << " " << y << "\n";
      //             std::cout << "Exception caught: " << e.what() << "\n";
      //             std::cout << "PRERHS: " << pre_rhs << "\n";
      //             std::cout << "PRELHS: " << pre_lhs << "\n\n";
      //             std::cout << "RHS: " << rhs << "\n";
      //             std::cout << "LHS: " << lhs << "\n\n";
      //             std::cout << "DEBUG: " << rhs(0,1) << "   " << rhs(1,0) << "\n\n";
      //             exit(0);
    }
    d += lhs;

    //          std::cout << "Update: " << lhs << "     " << d << "     " << sqrt(error_total) << "    " << (sqrt(lhs[2]*lhs[2]+lhs[5]*lhs[5])) << "\n";

    // Termination condition
    if (norm_2(lhs) < 0.01) 
      break;
  }
  //        std::cout << "----> " << d << "\n\n";
  offset(0) = d[2];
  offset(1) = d[5];
      
  if (norm_2(offset) > AFFINE_SUBPIXEL_MAX_TRANSLATION || 
      offset(0) != offset(0) ||  // Check to make sure the offset is not NaN...
      offset(1) != offset(1) ) { // ... ditto.
    return PixelDisparity<float>();
  } else {
    PixelDisparity<float> returnval = m_disparity_map(x,y);
    returnval.h() += offset(0);
    returnval.v() += offset(1);
    return returnval;
  }
}

void AffineSubpixelView::cache(BBox2i bbox) {

  //      std::cout << "bbox: " << bbox << "\n";
  m_left_bbox = bbox;
  m_left_bbox.min() -= Vector2i(m_kern_width/2+1,m_kern_height/2+1);
  m_left_bbox.max() += Vector2i(m_kern_width/2+1,m_kern_height/2+1);
  //      std::cout << "left bbox: " << m_left_bbox << "\n";

  ImageView<float> left_buf = crop( m_left_log_image, m_left_bbox );
  ImageView<float> x_deriv_buf = derivative_filter(crop( m_left_log_image, m_left_bbox ), 1, 0);
  ImageView<float> y_deriv_buf = derivative_filter(crop( m_left_log_image, m_left_bbox ), 0, 1);
  m_left_cached_log_image.reset(new CropView<ImageView<float> >( left_buf,
                                                                 BBox2i(-m_left_bbox.min().x(), -m_left_bbox.min().y(),
                                                                        m_left_image.cols(), m_left_image.rows()) ) );

  m_x_deriv.reset(new CropView<ImageView<float> >( x_deriv_buf,
                                                   BBox2i(-m_left_bbox.min().x(), -m_left_bbox.min().y(),
                                                          m_left_image.cols(), m_left_image.rows()) ) );

  m_y_deriv.reset(new CropView<ImageView<float> >( y_deriv_buf,
                                                   BBox2i(-m_left_bbox.min().x(), -m_left_bbox.min().y(),
                                                          m_left_image.cols(), m_left_image.rows()) ) );

  int num_good;
  BBox2 disp_range = disparity::get_disparity_range(crop(edge_extend(m_disparity_map,ZeroEdgeExtension()), m_left_bbox), num_good, false);
  //      std::cout << "disparity range: " << disp_range << "\n";

  m_right_bbox = bbox;
  m_right_bbox.min() -= Vector2i(m_kern_width/2+1,m_kern_height/2+1);
  m_right_bbox.min() += disp_range.min();
  m_right_bbox.max() += Vector2i(m_kern_width/2+1,m_kern_height/2+1);
  m_right_bbox.max() += disp_range.max();
  //      std::cout << "Right bbox: " << m_right_bbox << "\n";

  ImageView<float> right_buf = crop( m_right_log_image, m_right_bbox );
  m_right_cached_log_image.reset(new CropView<ImageView<float> >(right_buf, BBox2i(-m_right_bbox.min().x(), -m_right_bbox.min().y(),
                                                                                   m_right_image.cols(), m_right_image.rows()) ) );

}


