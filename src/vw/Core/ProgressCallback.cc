// __BEGIN_LICENSE__
// 
// Copyright (C) 2006 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration
// (NASA).  All Rights Reserved.
// 
// Copyright 2006 Carnegie Mellon University. All rights reserved.
// 
// This software is distributed under the NASA Open Source Agreement
// (NOSA), version 1.3.  The NOSA has been approved by the Open Source
// Initiative.  See the file COPYING at the top of the distribution
// directory tree for the complete NOSA document.
// 
// THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF ANY
// KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT
// LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL CONFORM TO
// SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR
// A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY WARRANTY THAT
// THE SUBJECT SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY THAT
// DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE SUBJECT SOFTWARE.
// 
// __END_LICENSE__

#include <vw/Core/ProgressCallback.h>

namespace {
  vw::ProgressCallback g_dummy_progress_callback_instance;
}

const vw::ProgressCallback &vw::ProgressCallback::dummy_instance() {
  return g_dummy_progress_callback_instance;
}

void vw::TerminalProgressCallback::report_progress(double progress) const {
  if (fabs(progress - m_last_reported_progress) > 0.01) {
    m_last_reported_progress = progress;
    int pi = static_cast<int>(progress * 60);
    vw_out(m_level) << "\r" << m_pre_progress_text << "[";
    for( int i=0; i<pi; ++i ) vw_out(m_level) << "*";
    for( int i=60; i>pi; --i ) vw_out(m_level) << ".";
    vw_out(m_level) << "] " << (int)(progress*100) << "%" << std::flush;
  }
}

void vw::TerminalProgressCallback::report_aborted(std::string why) const {
  vw_out(m_level) << " Aborted: " << why << std::endl;
}

void vw::TerminalProgressCallback::report_finished() const {
  vw_out(m_level) << "\r" << m_pre_progress_text << "[************************************************************] Complete!" << std::endl;
}
