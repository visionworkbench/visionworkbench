// __BEGIN_LICENSE__
// Copyright (C) 2006-2009 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


//template vw::TransformView<vw::InterpolationView<vw::EdgeExtensionView<ViewT, vw::ZeroEdgeExtension>, vw::BilinearInterpolation>, vw::cartography::GeoTransform> vw::cartography::geo_transform < ImageT > (const vw::ImageViewBase<ImageT>&, const vw::cartography::GeoReference&, const vw::cartography::GeoReference&);
//template vw::TransformView<vw::InterpolationView<vw::EdgeExtensionView<ViewT, vw::ZeroEdgeExtension>, vw::BilinearInterpolation>, vw::cartography::GeoTransform> vw::cartography::geo_transform < ImageT > (const vw::ImageViewBase<ImageT>&, const vw::cartography::GeoReference&, const vw::cartography::GeoReference&, vw::int32, vw::int32);
//template vw::TransformView<vw::InterpolationView<vw::EdgeExtensionView<ImageT, ExtensionT>, vw::BilinearInterpolation>, vw::cartography::GeoTransform> vw::cartography::geo_transform < ImageT,EdgeT > (const vw::ImageViewBase<ImageT>&, const vw::cartography::GeoReference&, const vw::cartography::GeoReference&, const EdgeT&);
//template vw::TransformView<vw::InterpolationView<vw::EdgeExtensionView<ImageT, ExtensionT>, vw::BilinearInterpolation>, vw::cartography::GeoTransform> vw::cartography::geo_transform < ImageT,EdgeT > (const vw::ImageViewBase<ImageT>&, const vw::cartography::GeoReference&, const vw::cartography::GeoReference&, vw::int32, vw::int32, const EdgeT&);
//template boost::disable_if<vw::IsScalar<InterpT>, vw::TransformView<vw::InterpolationView<vw::EdgeExtensionView<ImageT, ExtensionT>, InterpT>, vw::cartography::GeoTransform> >::type vw::cartography::geo_transform < ImageT,EdgeT,InterpT > (const vw::ImageViewBase<ImageT>&, const vw::cartography::GeoReference&, const vw::cartography::GeoReference&, const EdgeT&, const InterpT&);
//template vw::TransformView<vw::InterpolationView<vw::EdgeExtensionView<ImageT, ExtensionT>, InterpT>, vw::cartography::GeoTransform> vw::cartography::geo_transform < ImageT,EdgeT,InterpT > (const vw::ImageViewBase<ImageT>&, const vw::cartography::GeoReference&, const vw::cartography::GeoReference&, vw::int32, vw::int32, const EdgeT&, const InterpT&);
//template vw::math::Vector<DstElemT, 3> vw::cartography::lon_lat_radius_to_xyz < ElemT > (const vw::math::Vector<DstElemT, 3>&, bool);
//template vw::UnaryPerPixelView<ImageT, vw::cartography::LonLatRadToXYZFunctor> vw::cartography::lon_lat_radius_to_xyz < ImageT > (const vw::ImageViewBase<ImageT>&, bool);
//template vw::cartography::OrthoImageView<TerrainImageT, CameraImageT, InterpT, EdgeT> vw::cartography::orthoproject < TerrainImageT,CameraImageT,InterpT,EdgeT > (const vw::ImageViewBase<ImageT>&, const vw::cartography::GeoReference&, const vw::ImageViewBase<View2T>&, boost::shared_ptr<vw::camera::CameraModel>, const InterpT&, const EdgeT&);
//template int32 vw::cartography::output::kml::compute_resolution < TransformT > (const TransformT&, const vw::Vector2&);
//template int vw::cartography::output::tms::compute_resolution < TransformT > (const TransformT&, const vw::Vector2&);
//template vw::UnaryPerPixelView<ImageT, vw::cartography::ProjectPointFunctor> vw::cartography::project_point_image < ImageT > (const vw::ImageViewBase<ImageT>&, const vw::cartography::GeoReference&);
//template void vw::cartography::read_georeferenced_image < PixelT > (vw::ImageView<PixelT>&, vw::cartography::GeoReference&, const std::string&);
//template vw::UnaryPerPixelView<ImageT, vw::cartography::ReprojectPointFunctor> vw::cartography::reproject_point_image < ImageT > (const vw::ImageViewBase<ImageT>&, const vw::cartography::GeoReference&, const vw::cartography::GeoReference&);
//template void vw::cartography::write_georeferenced_image < ImageT > (const std::string&, const vw::ImageViewBase<ImageT>&, const vw::cartography::GeoReference&, const vw::ProgressCallback&);
//template vw::math::Vector<DstElemT, 3> vw::cartography::xyz_to_lon_lat_radius < ElemT > (const vw::math::Vector<DstElemT, 3>&, bool);
//template vw::UnaryPerPixelView<ImageT, vw::cartography::XYZtoLonLatRadFunctor> vw::cartography::xyz_to_lon_lat_radius < ImageT > (const vw::ImageViewBase<ImageT>&, bool);
