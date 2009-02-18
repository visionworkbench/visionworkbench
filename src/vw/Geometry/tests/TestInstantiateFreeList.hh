// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


//template bool vw::geometry::operator!= < ShapeT,RealT,DimN > (const vw::geometry::ShapeBase<ShapeT, RealT, DimN>&, const vw::geometry::ShapeBase<ShapeT, RealT, DimN>&);
//template bool vw::geometry::operator!= < SphereT1,RealT1,DimN1,SphereT2,RealT2,DimN2 > (const vw::geometry::SphereBase<SphereT, RealT, DimN>&, const vw::geometry::SphereBase<SphereT2, RealT2, DimN2>&);
//template vw::geometry::Box<RealT, DimN> vw::geometry::operator* < RealT,DimN,ScalarT > (ScalarT, const vw::geometry::Box<RealT, DimN>&);
//template vw::geometry::Box<RealT, DimN> vw::geometry::operator* < RealT,DimN,ScalarT > (const vw::geometry::Box<RealT, DimN>&, ScalarT);
//template ShapeT vw::geometry::operator* < ShapeT,RealT,DimN,ScalarT > (ScalarT, const vw::geometry::ShapeBase<ShapeT, RealT, DimN>&);
//template ShapeT vw::geometry::operator* < ShapeT,RealT,DimN,ScalarT > (const vw::geometry::ShapeBase<ShapeT, RealT, DimN>&, ScalarT);
//template SphereT vw::geometry::operator* < SphereT,RealT,DimN,ScalarT > (ScalarT, const vw::geometry::SphereBase<SphereT, RealT, DimN>&);
//template SphereT vw::geometry::operator* < SphereT,RealT,DimN,ScalarT > (const vw::geometry::SphereBase<SphereT, RealT, DimN>&, ScalarT);
//template vw::geometry::Box<RealT, DimN> vw::geometry::operator+ < RealT,DimN,VectorT > (const vw::geometry::Box<RealT, DimN>&, const vw::math::VectorBase<VectorT>&);
//template vw::geometry::Box<RealT, DimN> vw::geometry::operator+ < RealT,DimN,VectorT > (const vw::math::VectorBase<VectorT>&, const vw::geometry::Box<RealT, DimN>&);
//template ShapeT vw::geometry::operator+ < ShapeT,RealT,DimN,VectorT > (const vw::geometry::ShapeBase<ShapeT, RealT, DimN>&, const vw::math::VectorBase<VectorT>&);
//template ShapeT vw::geometry::operator+ < ShapeT,RealT,DimN,VectorT > (const vw::math::VectorBase<VectorT>&, const vw::geometry::ShapeBase<ShapeT, RealT, DimN>&);
//template SphereT vw::geometry::operator+ < SphereT,RealT,DimN,VectorT > (const vw::geometry::SphereBase<SphereT, RealT, DimN>&, const vw::math::VectorBase<VectorT>&);
//template SphereT vw::geometry::operator+ < SphereT,RealT,DimN,VectorT > (const vw::math::VectorBase<VectorT>&, const vw::geometry::SphereBase<SphereT, RealT, DimN>&);
//template vw::geometry::Box<RealT, DimN> vw::geometry::operator- < RealT,DimN,VectorT > (const vw::geometry::Box<RealT, DimN>&, const vw::math::VectorBase<VectorT>&);
//template ShapeT vw::geometry::operator- < ShapeT,RealT,DimN,VectorT > (const vw::geometry::ShapeBase<ShapeT, RealT, DimN>&, const vw::math::VectorBase<VectorT>&);
//template SphereT vw::geometry::operator- < SphereT,RealT,DimN,VectorT > (const vw::geometry::SphereBase<SphereT, RealT, DimN>&, const vw::math::VectorBase<VectorT>&);
//template vw::geometry::Box<RealT, DimN> vw::geometry::operator/ < RealT,DimN,ScalarT > (const vw::geometry::Box<RealT, DimN>&, ScalarT);
//template ShapeT vw::geometry::operator/ < ShapeT,RealT,DimN,ScalarT > (const vw::geometry::ShapeBase<ShapeT, RealT, DimN>&, ScalarT);
//template SphereT vw::geometry::operator/ < SphereT,RealT,DimN,ScalarT > (const vw::geometry::SphereBase<SphereT, RealT, DimN>&, ScalarT);
//template bool vw::geometry::operator== < ShapeT,RealT,DimN > (const vw::geometry::ShapeBase<ShapeT, RealT, DimN>&, const vw::geometry::ShapeBase<ShapeT, RealT, DimN>&);
//template bool vw::geometry::operator== < SphereT1,RealT1,DimN1,SphereT2,RealT2,DimN2 > (const vw::geometry::SphereBase<SphereT, RealT, DimN>&, const vw::geometry::SphereBase<SphereT2, RealT2, DimN2>&);
//template void vw::geometry::read_box < RealT,DimN > (const std::string&, vw::geometry::Box<RealT, DimN>&, bool);
//template void vw::geometry::write_box < RealT,DimN > (const std::string&, const vw::geometry::Box<RealT, DimN>&, bool);
//template void vw::math::read_point_list < ContainerT > (const std::string&, std::vector<CharT, std::allocator<_CharT> >&, bool);
//template void vw::math::read_point_list < ContainerT > (std::istream&, std::vector<CharT, std::allocator<_CharT> >&, bool);
//template void vw::math::write_point_list < ContainerT > (const std::string&, const std::vector<CharT, std::allocator<_CharT> >&, bool);
//template void vw::math::write_point_list < ContainerT > (std::ostream&, const std::vector<CharT, std::allocator<_CharT> >&, bool);
