// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <cxxtest/TestSuite.h>
#include <string>
#include <vw/Core.h>
#include <vw/Image/PixelTypes.h>
#include <vw/Stereo/DisparityMap.h>
#include <boost/shared_ptr.hpp>

using namespace vw;

const bool Arg1IsCompound = false;
const bool Arg2IsCompound = false;
const bool CompoundB = true;
const int ChannelsN = 3;

typedef int T;
typedef int T1;
typedef int T2;

typedef int DstT;
typedef int SrcT;
typedef int Val1T;
typedef int Val2T;
typedef int Val3T;
typedef int ValT;

typedef short Arg1T[];
typedef int Arg2T[];
typedef int Args;
typedef int ArgT[];
typedef int ArgsT;
typedef int ChannelT;
typedef char CharT;

typedef std::char_traits<char> traits;

struct GeneratorT {
  typedef int value_type;
  int size() const {}
  boost::shared_ptr<int> generate() {}
};

struct TaskT {
  void operator()() {}
};

struct FuncT : ReturnFixedType<int> {
  template <class Args>
  int operator() (int a) {return a;}
  int operator() (int a, int b) {return a+b;}
};

typedef int FuncArg1T;
typedef int FuncArg2T;
typedef int FuncArgT;

struct ResultT
{
  template <typename T>
  ResultT(T a=0, T b=0, T c=0, T d=0) {}
};

#include "TestInstantiateRecordList.hh"

class TestInstantiateCoreRecord : public CxxTest::TestSuite
{
  public: void test_inst() {}
};
