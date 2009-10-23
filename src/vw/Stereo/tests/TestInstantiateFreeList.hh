// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


//template vw::PixelDisparity<float> vw::stereo::compute_disparity < ChannelT > (vw::ImageView<PixelT>&, vw::ImageView<PixelT>&, int, int, int, int, int, int, int, int);
//template double vw::stereo::compute_soad < ChannelT > (ChannelT*, ChannelT*, int, int, int, int, int, int, int, int);
//template double vw::stereo::compute_soad < ViewT > (const vw::ImageViewBase<ImageT>&, const vw::ImageViewBase<ImageT>&, int, int, int, int, int, int, const vw::BBox2i&, const vw::BBox2i&);
//template vw::UnaryPerPixelAccessorView<vw::EdgeExtensionView<vw::UnaryPerPixelAccessorView<vw::EdgeExtensionView<ViewT, vw::ZeroEdgeExtension>, vw::stereo::disparity::RemoveOutliersFunc>, vw::ZeroEdgeExtension>, vw::stereo::disparity::RemoveOutliersFunc> vw::stereo::disparity::clean_up < ViewT > (const vw::ImageViewBase<ImageT>&, int, int, double, double);
//template vw::ImageView<bool> vw::stereo::disparity::generate_mask < ViewT > (const vw::ImageViewBase<ImageT>&, vw::int32);
//template BBox2 vw::stereo::disparity::get_disparity_range < ViewT > (const vw::ImageViewBase<ImageT>&, bool);
//template BBox2 vw::stereo::disparity::get_disparity_range < ViewT > (const vw::ImageViewBase<ImageT>&, int&, bool, const vw::ProgressCallback&);
//template vw::UnaryPerPixelView<ViewT, vw::stereo::disparity::LessThanThresholdFunc> vw::stereo::disparity::less_than_threshold < ViewT > (const vw::ImageViewBase<ImageT>&, double);
//template vw::BinaryPerPixelView<ViewT, vw::PixelIndex3View, vw::stereo::disparity::DisparityMaskFunc<MaskViewT> > vw::stereo::disparity::mask < ViewT,MaskViewT > (const vw::ImageViewBase<ImageT>&, const vw::ImageViewBase<View2T>&, const vw::ImageViewBase<View2T>&);
//template void vw::stereo::disparity::mask_black_pixels < ViewT > (const vw::ImageViewBase<ImageT>&, vw::ImageView<unsigned char>&);
//template vw::UnaryPerPixelView<ViewT, vw::stereo::disparity::MissingPixelImageFunc> vw::stereo::disparity::missing_pixel_image < ViewT > (vw::ImageViewBase<ImageT>&);
//template vw::BinaryPerPixelView<ViewT, vw::PixelIndex3View, vw::stereo::disparity::BorderPixelsFunc> vw::stereo::disparity::remove_border_pixels < ViewT > (const vw::ImageViewBase<ImageT>&, int);
//template vw::BinaryPerPixelView<ViewT, vw::PixelIndex3View, vw::stereo::disparity::InvalidPixelsFunc> vw::stereo::disparity::remove_invalid_pixels < ViewT > (vw::ImageViewBase<ImageT>&, int, int);
//template vw::UnaryPerPixelAccessorView<vw::EdgeExtensionView<ViewT, vw::ZeroEdgeExtension>, vw::stereo::disparity::RemoveOutliersFunc> vw::stereo::disparity::remove_outliers < ViewT > (const vw::ImageViewBase<ImageT>&, int, int, double, double);
//template vw::UnaryPerPixelAccessorView<vw::EdgeExtensionView<ViewT, vw::ZeroEdgeExtension>, vw::stereo::disparity::StdDevImageFunc> vw::stereo::disparity::std_dev_image < ViewT > (const vw::ImageViewBase<ImageT>&, int, int);
//template vw::UnaryPerPixelAccessorView<vw::EdgeExtensionView<ImageT, ExtensionT>, vw::stereo::disparity::StdDevImageFunc> vw::stereo::disparity::std_dev_image < ViewT,EdgeT > (const vw::ImageViewBase<ImageT>&, int, int, EdgeT);
//template vw::BinaryPerPixelView<ViewT, vw::PixelIndex3View, vw::stereo::disparity::TransformDisparitiesFunc<TransformT> > vw::stereo::disparity::transform_disparities < ViewT,TransformT > (const vw::ImageViewBase<ImageT>&, const TransformT&);
//template void vw::stereo::subpixel_correlation_affine_2d < ChannelT > (vw::ImageView<vw::PixelDisparity<float> >&, const vw::ImageView<PixelT>&, const vw::ImageView<PixelT>&, int, int, bool, bool, bool);
//template void vw::stereo::subpixel_correlation_affine_2d_EM < ChannelT > (vw::ImageView<vw::PixelDisparity<float> >&, const vw::ImageView<PixelT>&, const vw::ImageView<PixelT>&, int, int, bool, bool, bool);
//template void vw::stereo::subpixel_correlation_affine_2d_bayesian < ChannelT > (vw::ImageView<vw::PixelDisparity<float> >&, const vw::ImageView<PixelT>&, const vw::ImageView<PixelT>&, int, int, bool, bool, bool);
//template void vw::stereo::subpixel_correlation_parabola < ChannelT > (vw::ImageView<vw::PixelDisparity<float> >&, const vw::ImageView<PixelT>&, const vw::ImageView<PixelT>&, int, int, bool, bool, bool);
