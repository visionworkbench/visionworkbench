#ifndef CHIPPER_H
#define CHIPPER_H

/******************************************************************************
 * $Id$
 *
 * Project:  PDAL - http://www.pdal.io/index.html - A BSD library for point cloud data.
 * Purpose:  PDAL chipper class
 * Author:   Howard Butler, hobu.inc@gmail.com
 *
 ******************************************************************************
 * Copyright (c) 2011, Howard Butler
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided
 *       with the distribution.
 *     * Neither the name of the Martin Isenburg or Iowa Department
 *       of Natural Resources nor the names of its contributors may be
 *       used to endorse or promote products derived from this software
 *       without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 ****************************************************************************/

#include <vw/Math/Vector.h>
#include <vw/Image/ImageView.h>
#include <vw/Cartography/GeoReference.h>


/**
The objective is to split the region into non-overlapping blocks, each
containing approximately the same number of points, as specified by the
user.  We'd also like the blocks closer to square than not.

**/

namespace pdal
{

//class Stage;

namespace filters

{

typedef boost::uint64_t PointId;
typedef boost::uint64_t point_count_t;
typedef std::vector<vw::Vector3> PointBuffer;

class Chipper;

enum Direction
{
    DIR_X,
    DIR_Y,
    DIR_NONE
};


class VW_API ChipPtRef
{
    friend class ChipRefList;
    friend class Chipper;

private:
    double m_pos;
    boost::uint32_t m_ptindex;
    boost::uint32_t m_oindex;

public:
    bool operator < (const ChipPtRef& pt) const
    {
        return m_pos < pt.m_pos;
    }
};

class VW_API ChipRefList
{
    friend class Chipper;

private:
    std::vector<ChipPtRef> m_vec;
    Direction m_dir;

    ChipRefList(Direction dir = DIR_NONE) : m_dir(dir)
    {}
    std::vector<ChipPtRef>::size_type size() const
    {
        return m_vec.size();
    }
    void reserve(std::vector<ChipPtRef>::size_type n)
    {
        m_vec.reserve(n);
    }
    void resize(std::vector<ChipPtRef>::size_type n)
    {
        m_vec.resize(n);
    }
    void push_back(const ChipPtRef& ref)
    {
        m_vec.push_back(ref);
    }
    std::vector<ChipPtRef>::iterator begin()
    {
        return m_vec.begin();
    }
    std::vector<ChipPtRef>::iterator end()
    {
        return m_vec.end();
    }
    ChipPtRef& operator[](boost::uint32_t pos)
    {
        return m_vec[pos];
    }
    std::string Dir()
    {
        if (m_dir == DIR_X)
            return "X";
        else if (m_dir == DIR_Y)
            return "Y";
        else
            return "NONE";
    }
};

class VW_API Chipper
{
public:

  Chipper(PointBuffer& buffer, int blockSize,
          bool have_georef, vw::cartography::GeoReference const& georef,
          int num_out_cols, int num_out_rows,
          vw::ImageView<vw::Vector3> & outImg);
  
private:
    void load(PointBuffer const& buffer, ChipRefList& xvec,
        ChipRefList& yvec, ChipRefList& spare);
    void partition(point_count_t size);
    void decideSplit(ChipRefList& v1, ChipRefList& v2,
        ChipRefList& spare, PointId left, PointId right);
    void split(ChipRefList& wide, ChipRefList& narrow,
        ChipRefList& spare, PointId left, PointId right);
    void finalSplit(ChipRefList& wide, ChipRefList& narrow,
        PointId pleft, PointId pcenter);
    void emit(ChipRefList& wide, PointId widemin, PointId widemax,
        ChipRefList& narrow, PointId narrowmin, PointId narrowmax);

    PointBuffer const& m_inbuf;
    int m_blockSize;
    bool m_have_georef;
    vw::cartography::GeoReference m_georef;
    int m_numMaxPtsInChip;
    int m_currChip;
    std::vector<PointId> m_partitions;
    ChipRefList m_xvec;
    ChipRefList m_yvec;
    ChipRefList m_spare;
    vw::ImageView<vw::Vector3> & m_outImg;
    
    Chipper& operator=(const Chipper&); // not implemented
    Chipper(const Chipper&); // not implemented
};

} // namespace filters
} // namespace pdal

#endif
