// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


//template void vw::camera::bsolve < T > (vw::math::Vector<ElemT, 0>&, vw::math::Matrix<ElemT, 0, 0>&);
//template vw::TransformView<vw::InterpolationView<vw::EdgeExtensionView<ViewT, vw::ZeroEdgeExtension>, vw::BilinearInterpolation>, vw::camera::CameraTransform<SrcCameraT, DstCameraT> > vw::camera::camera_transform < ImageT,SrcCameraT,DstCameraT > (const vw::ImageViewBase<ImageT>&, const SrcCameraT&, const DstCameraT&);
//template vw::TransformView<vw::InterpolationView<vw::EdgeExtensionView<ImageT, EdgeT>, vw::BilinearInterpolation>, vw::camera::CameraTransform<SrcCameraT, DstCameraT> > vw::camera::camera_transform < ImageT,SrcCameraT,DstCameraT,EdgeT > (const vw::ImageViewBase<ImageT>&, const SrcCameraT&, const DstCameraT&, const EdgeT&);
//template vw::TransformView<vw::InterpolationView<vw::EdgeExtensionView<ImageT, EdgeT>, InterpT>, vw::camera::CameraTransform<SrcCameraT, DstCameraT> > vw::camera::camera_transform < ImageT,SrcCameraT,DstCameraT,EdgeT,InterpT > (const vw::ImageViewBase<ImageT>&, const SrcCameraT&, const DstCameraT&, const EdgeT&, const InterpT&);
//template unsigned int vw::camera::chol_inverse < T > (vw::math::Matrix<ElemT, 0, 0>&);
//template void vw::camera::cholesky < T > (vw::math::Matrix<ElemT, 0, 0>&);
//template void vw::camera::fsolve < T > (vw::math::Vector<ElemT, 0>&, vw::math::Matrix<ElemT, 0, 0>&);
//template vw::ImageView<vw::PixelRGB<vw::CompoundChannelType<ImplT::pixel_type>::type> > vw::camera::inverse_bayer_filter < ViewT > (const vw::ImageViewBase<ImageT>&);
//template void vw::camera::ldl_decomposition < ElemT > (vw::camera::SparseSkylineMatrix<ElemT>&);
//template vw::TransformView<vw::InterpolationView<vw::EdgeExtensionView<ViewT, vw::ZeroEdgeExtension>, vw::BilinearInterpolation>, vw::camera::CameraTransform<SrcCameraT, SrcCameraT::linearized_type> > vw::camera::linearize_camera_transform < ImageT,SrcCameraT > (const vw::ImageViewBase<ImageT>&, const SrcCameraT&);
//template vw::TransformView<vw::InterpolationView<vw::EdgeExtensionView<ViewT, EdgeT>, vw::BilinearInterpolation>, vw::camera::CameraTransform<SrcCameraT, SrcCameraT::linearized_type> > vw::camera::linearize_camera_transform < ImageT,SrcCameraT,EdgeT > (const vw::ImageViewBase<ImageT>&, const SrcCameraT&, const EdgeT&);
//template vw::TransformView<vw::InterpolationView<vw::EdgeExtensionView<ViewT, EdgeT>, InterpT>, vw::camera::CameraTransform<SrcCameraT, SrcCameraT::linearized_type> > vw::camera::linearize_camera_transform < ImageT,SrcCameraT,EdgeT,InterpT > (const vw::ImageViewBase<ImageT>&, const SrcCameraT&, const EdgeT&, const InterpT&);
//template unsigned int vw::camera::mod_cholesky < T > (vw::math::Matrix<ElemT, 0, 0>&);
//template unsigned int vw::camera::solve < T > (vw::math::Vector<ElemT, 0>&, vw::math::Matrix<ElemT, 0, 0>&);
//template vw::math::Vector<double, 0> vw::camera::sparse_solve < ElemT,VectorT > (vw::camera::SparseSkylineMatrix<ElemT>&, const VectorT&);
