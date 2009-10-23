// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


//template boost::enable_if<boost::is_integral<NumberT>, float>::type vw::ip::HHaarWavelet < ViewT,NumberT > (const vw::ImageViewBase<ImageT>&, const NumberT&, const NumberT&, const float&);
//template boost::enable_if<boost::is_floating_point<NumberT>, float>::type vw::ip::HHaarWavelet < ViewT,NumberT > (const vw::ImageViewBase<ImageT>&, const NumberT&, const NumberT&, const float&);
//template vw::ImageView<double> vw::ip::IntegralImage < ViewT > (const vw::ImageViewBase<ImageT>&);
//template boost::enable_if<boost::is_integral<NumberT>, float>::type vw::ip::VHaarWavelet < ViewT,NumberT > (const vw::ImageViewBase<ImageT>&, const NumberT&, const NumberT&, const float&);
//template boost::enable_if<boost::is_floating_point<NumberT>, float>::type vw::ip::VHaarWavelet < ViewT,NumberT > (const vw::ImageViewBase<ImageT>&, const NumberT&, const NumberT&, const float&);
//template InterestPointList vw::ip::crop < RealT > (const vw::ip::InterestPointList&, const vw::math::BBox<RealT, 2>&);
//template InterestPointList vw::ip::detect_interest_points < ViewT,DetectorT > (const ViewT&, DetectorT&, int);
//template int vw::ip::filter_1d < T1,T2 > (std::vector<CharT, std::allocator<_CharT> >&, std::vector<DomainT, std::allocator<_T2> >&, bool);
//template int vw::ip::find_peaks < DataT > (vw::ip::InterestPointList&, const DataT&);
//template int vw::ip::find_peaks < DataT,ViewT > (vw::ip::InterestPointList&, const std::vector<CharT, std::allocator<_CharT> >&, const vw::ip::ImageOctave<ViewT>&);
//template void vw::ip::find_weighted_histogram_mode < T > (const std::vector<CharT, std::allocator<_CharT> >&, std::vector<int, std::allocator<int> >&);
//template bool vw::ip::fit_peak < DataT > (const std::vector<CharT, std::allocator<_CharT> >&, vw::ip::InterestPoint&, const vw::ip::ImageOctave<DataT::source_type>&, float*, float*, float*);
//template bool vw::ip::fit_peak < ViewT > (const vw::ImageViewBase<ImageT>&, vw::ip::InterestPoint&, float*, float*);
//template int vw::ip::fit_peak_1D < ElmtT > (const vw::math::Vector<DstElemT, 3>&, ElmtT&, ElmtT*);
//template float vw::ip::get_orientation < OriT,MagT > (const vw::ImageViewBase<ImageT>&, const vw::ImageViewBase<View2T>&, float, float, float);
//template bool vw::ip::is_extremum < DataT > (const std::vector<CharT, std::allocator<_CharT> >&, int, int, int, int);
//template bool vw::ip::is_extremum < ViewT > (const vw::ImageViewBase<ImageT>&, int, int, int);
//template bool vw::ip::is_local_max < DataT > (const std::vector<CharT, std::allocator<_CharT> >&, int, int, int);
//template bool vw::ip::is_local_max < ViewT > (const vw::ImageViewBase<ImageT>&, int, int);
//template bool vw::ip::is_local_min < DataT > (const std::vector<CharT, std::allocator<_CharT> >&, int, int, int);
//template bool vw::ip::is_local_min < ViewT > (const vw::ImageViewBase<ImageT>&, int, int);
//template bool vw::ip::is_local_minmax < DataT > (const std::vector<CharT, std::allocator<_CharT> >&, int, int, int);
//template bool vw::ip::is_local_minmax < ViewT > (const vw::ImageViewBase<ImageT>&, int, int);
//template int vw::ip::make_gaussian_kernel_2d < KernelT > (KernelT&, float, int);
//template void vw::ip::non_max_suppression < T > (std::vector<CharT, std::allocator<_CharT> >&, bool);
//template void vw::ip::orientation_histogram < ViewT1,ViewT2 > (const vw::ImageViewBase<ImageT>&, const vw::ImageViewBase<View2T>&, std::vector<float, std::allocator<float> >&, int, int, float, unsigned int);
//template int vw::ip::smooth_weighted_histogram < T > (std::vector<CharT, std::allocator<_CharT> >&, float);
//template void vw::ip::weighted_histogram < ViewT1,ViewT2 > (const vw::ImageViewBase<ImageT>&, const vw::ImageViewBase<View2T>&, std::vector<float, std::allocator<float> >&, float, float, unsigned int);
//template void vw::ip::weighted_magnitude < ViewT1,ViewT2 > (vw::ImageViewBase<ImageT>&, const vw::ImageViewBase<View2T>&, int, int, float, int);
//template void vw::read_matrix < T > (vw::math::Matrix<ElemT, 0, 0>&, const std::string&);
//template void vw::read_vector < T > (vw::math::Vector<ElemT, 0>&, const std::string&);
//template void vw::write_matrix < T > (const std::string&, vw::math::Matrix<ElemT, 0, 0>&);
//template void vw::write_vector < T > (const std::string&, vw::math::Vector<ElemT, 0>&);
