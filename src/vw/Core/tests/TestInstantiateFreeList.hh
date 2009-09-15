// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


template T vw::_numeric < T > (T);
template vw::CompoundResult<FuncT, Arg1T, Arg2T>::type vw::compound_apply < FuncT,Arg1T,Arg2T > (const FuncT&, const Arg1T&, const Arg2T&);
template vw::CompoundResult<FuncT, ArgT, void>::type vw::compound_apply < FuncT,ArgT > (const FuncT&, const ArgT&);
//template ResultT vw::compound_select_channel < ResultT,PixelT > (PixelT&, boost::enable_if<boost::mpl::and_<boost::mpl::not_<vw::IsCompound<PixelT> >::type, boost::is_reference<T>::type, mpl_::bool_<true>, mpl_::bool_<true>, mpl_::bool_<true> >, int>::type);
//template ResultT vw::compound_select_channel < ResultT,PixelT > (PixelT&, boost::enable_if<boost::mpl::and_<vw::IsCompound<PixelT>, boost::is_reference<T>::type, mpl_::bool_<true>, mpl_::bool_<true>, mpl_::bool_<true> >, int>::type);
//template ResultT vw::compound_select_channel < ResultT,PixelT > (const PixelT&, boost::enable_if<boost::mpl::and_<boost::mpl::not_<vw::IsCompound<PixelT> >::type, boost::mpl::not_<boost::is_reference<T>::type>::type, mpl_::bool_<true>, mpl_::bool_<true>, mpl_::bool_<true> >, int>::type);
template ResultT vw::compound_select_channel < ResultT,PixelT > (const PixelT&, boost::enable_if<boost::mpl::and_<vw::IsCompound<PixelT>, boost::mpl::not_<boost::is_reference<T>::type>::type, mpl_::bool_<true>, mpl_::bool_<true>, mpl_::bool_<true> >, int>::type);
template boost::enable_if<vw::IsScalarOrCompound<T>, double>::type vw::mean_channel_value < T > (const T&);
template std::string vw::stringify < T > (const T&);
