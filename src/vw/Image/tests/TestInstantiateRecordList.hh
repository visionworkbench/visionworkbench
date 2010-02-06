// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


template class  vw::BinaryCompoundFunctor<FuncT, vw::PixelMask<ChildPixel1T>, vw::PixelMask<ChildPixel2T> >;
template class  vw::BinaryInPlaceCompoundFunctor<FuncT, vw::PixelMask<ChildPixel1T>, vw::PixelMask<ChildPixel2T> >;
template struct vw::ChannelRangeHelper<T, false>;
template struct vw::CompoundChannelCast<vw::PixelGray<ChT>, NewChT>;
template struct vw::CompoundChannelCast<vw::PixelGray<ChT>, const NewChT>;
template struct vw::CompoundChannelCast<vw::PixelGrayA<ChT>, NewChT>;
template struct vw::CompoundChannelCast<vw::PixelGrayA<ChT>, const NewChT>;
template struct vw::CompoundChannelCast<vw::PixelHSV<ChT>, NewChT>;
template struct vw::CompoundChannelCast<vw::PixelHSV<ChT>, const NewChT>;
template struct vw::CompoundChannelCast<vw::PixelLab<ChT>, NewChT>;
template struct vw::CompoundChannelCast<vw::PixelLab<ChT>, const NewChT>;
template struct vw::CompoundChannelCast<vw::PixelLuv<ChT>, NewChT>;
template struct vw::CompoundChannelCast<vw::PixelLuv<ChT>, const NewChT>;
template struct vw::CompoundChannelCast<vw::PixelMask<ChildT>, NewChT>;
template struct vw::CompoundChannelCast<vw::PixelMask<ChildT>, const NewChT>;
template struct vw::CompoundChannelCast<vw::PixelRGB<ChT>, NewChT>;
template struct vw::CompoundChannelCast<vw::PixelRGB<ChT>, const NewChT>;
template struct vw::CompoundChannelCast<vw::PixelRGBA<ChT>, NewChT>;
template struct vw::CompoundChannelCast<vw::PixelRGBA<ChT>, const NewChT>;
template struct vw::CompoundChannelCast<vw::PixelXYZ<ChT>, NewChT>;
template struct vw::CompoundChannelCast<vw::PixelXYZ<ChT>, const NewChT>;
template struct vw::CompoundChannelCast<vw::math::Vector<OldChT, SizeN>, NewChT>;
template struct vw::CompoundChannelCast<vw::math::Vector<OldChT, SizeN>, const NewChT>;
template class  vw::UnaryCompoundFunctor<FuncT, vw::PixelMask<ChildPixel1T> >;
template class  vw::UnaryInPlaceCompoundFunctor<FuncT, vw::PixelMask<ChildPixel1T> >;
