template void vw::read_image < PixelT > (vw::ImageView<PixelT>&, const std::string&);
//template void vw::write_image < ElemT > (const std::string&, const std::vector<CharT, std::allocator<_CharT> >&);
template void vw::write_image < ImageT > (const std::string&, const vw::ImageViewBase<ImageT>&, const vw::ProgressCallback&);
