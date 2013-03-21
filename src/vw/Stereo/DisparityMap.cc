#include <vw/Stereo/DisparityMap.h>

namespace vw { namespace stereo {

  ImageView< PixelMask<Vector2i> >
  transform_disparities(BBox2i subregion,
                        Matrix<double> const& T,
                        ImageView< PixelMask<Vector2i> > const& disparity
                        ){

    VW_ASSERT(subregion.width() == disparity.cols() &&
              subregion.height() == disparity.rows(),
              ArgumentErr() << "transform_disparities: "
              << "The sizes of subregion and disparity don't match.\n");

    typedef PixelMask<Vector2i> IntDispT;
    ImageView<IntDispT> out_disp(disparity.cols(), disparity.rows());
    HomographyTransform H(T);

    for (int x = 0; x < disparity.cols(); x++){
      for (int y = 0; y < disparity.rows(); y++){

        IntDispT disp = disparity(x, y);
        if (!is_valid(disp)) continue; // output stays invalid

        Vector2 beg = subregion.min() + Vector2(x, y);
        Vector2 end = beg + disp.child();
        Vector2 end_trans = H.forward(end);

        out_disp(x, y) = round( end_trans - beg );
        out_disp(x, y).validate(); // unnecessary but clearer

      }
    }

    return out_disp;
  }

}}    // namespace vw::stereo
