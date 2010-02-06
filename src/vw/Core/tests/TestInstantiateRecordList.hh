// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


template struct vw::CompoundChannelCast<const T, ChannelT>;
template struct vw::CompoundChannelType<const T>;
template struct vw::CompoundIsCompatible<T1, const T2>;
template struct vw::CompoundNumChannels<const T>;
template struct vw::CompoundResult<FuncT, ArgT, void>;
template struct vw::IsCompound<const T>;
template struct vw::IsScalar<const T>;
template struct vw::PromoteTypeSpecializationHelper<Arg1T, Arg2T, false, true>;
template struct vw::PromoteTypeSpecializationHelper<Arg1T, Arg2T, true, false>;
template struct vw::PromoteTypeSpecializationHelper<Arg1T, Arg2T, true, true>;
template struct vw::TypeDeductionHelper<ArgT, ArgT>;
template struct vw::TypeDeductionHelper<T, double>;
template struct vw::TypeDeductionHelper<T, float>;
template struct vw::TypeDeductionHelper<T, int>;
template struct vw::TypeDeductionHelper<T, long double>;
template struct vw::TypeDeductionHelper<T, long int>;
template struct vw::TypeDeductionHelper<T, long unsigned int>;
template struct vw::TypeDeductionHelper<T, short int>;
template struct vw::TypeDeductionHelper<T, short unsigned int>;
template struct vw::TypeDeductionHelper<T, signed char>;
template struct vw::TypeDeductionHelper<T, unsigned char>;
template struct vw::TypeDeductionHelper<T, unsigned int>;
template struct vw::TypeDeductionHelper<double, T>;
template struct vw::TypeDeductionHelper<float, T>;
template struct vw::TypeDeductionHelper<int, T>;
template struct vw::TypeDeductionHelper<long double, T>;
template struct vw::TypeDeductionHelper<long int, T>;
template struct vw::TypeDeductionHelper<long unsigned int, T>;
template struct vw::TypeDeductionHelper<short int, T>;
template struct vw::TypeDeductionHelper<short unsigned int, T>;
template struct vw::TypeDeductionHelper<signed char, T>;
template struct vw::TypeDeductionHelper<unsigned char, T>;
template struct vw::TypeDeductionHelper<unsigned int, T>;
