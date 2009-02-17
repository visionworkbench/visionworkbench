// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


template void vw::read_image < PixelT > (vw::ImageView<PixelT>&, const std::string&);
//template void vw::write_image < ElemT > (const std::string&, const std::vector<CharT, std::allocator<_CharT> >&);
template void vw::write_image < ImageT > (const std::string&, const vw::ImageViewBase<ImageT>&, const vw::ProgressCallback&);
