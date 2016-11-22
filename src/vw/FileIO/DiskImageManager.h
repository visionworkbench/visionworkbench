// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NASA Vision Workbench is licensed under the Apache License,
//  Version 2.0 (the "License"); you may not use this file except in
//  compliance with the License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// A thread-safe way for keeping at most a given number of open handles
/// to a set of georeferenced images. When a handle for a new image is
/// is requested, at a certain spatial location, and we are at the
/// maximum number of open handles, close the handle to the image that
/// is spatially furthest from the desired location and open the
/// current one.

#ifndef __VW_FILEIO_DISKIMAGEMANAGER_H__
#define __VW_FILEIO_DISKIMAGEMANAGER_H__

#include <vw/FileIO/DiskImageView.h>

#include <vw/Core/Exception.h>
#include <vw/Core/Thread.h>
#include <vw/Core/FundamentalTypes.h>

#include <typeinfo>
#include <sstream>
#include <string>

#include <boost/smart_ptr/shared_ptr.hpp>

namespace vw {

  template<class PixelT>
  class DiskImageManager {
  private:
    typedef  boost::shared_ptr< DiskImageView<PixelT> > ImgPtr;
    std::vector<std::string> m_fileNames;
    std::vector<ImgPtr>      m_fileHandles;
    std::vector<BBox2>       m_imgBoxes;
    Mutex                    m_mutex;

    class DiskImageManagerLineBase;
    template <class GeneratorT> class DiskImageManagerLine;
    std::map<int, int> m_used; // open, and actively being used, can't be closed
    int m_max_num_open;
    bool m_handles_almost_used_up;
  public:
    // Do some forward declarations so we can put the public interface first
    DiskImageManager(): m_max_num_open(0), m_handles_almost_used_up(false){}
    size_t size(){
      return m_fileNames.size();
    }

    // We can't open more handles to files. Free up a few of them.
    // Later on, when we get some breathing room, we'll manage opening
    // and closing properly.
    void freeup_handles_not_thread_safe(){
      vw_out() << "Too many files were open. Will start selectively closing them.\n";
      m_handles_almost_used_up = true;
      m_max_num_open = m_fileHandles.size();
      m_max_num_open -= 100;
      if (m_max_num_open < 1) m_max_num_open = 1; // have at least one open handle
      for (int c = m_max_num_open; c < (int)m_fileHandles.size(); c++) {
	m_fileHandles[c] = ImgPtr();
      }
    }
      
    // This function is not thread-safe, by design. It won't be used
    // from multiple threads, and there is no need for a lock. Add a
    // file to be managed. Add a "box" in the same coordinate system
    // as later used in get_handle() to show where this file is
    // located in space.
    void add_file_handle_not_thread_safe(std::string const& fileName, BBox2 const& img_box){

      m_fileNames.push_back(fileName);
      m_imgBoxes.push_back(img_box);
      
      if (!m_handles_almost_used_up) {
	try {
          // Try to open a handle
	  vw_out() << "Loading: " << fileName << std::endl;
	  m_fileHandles.push_back( ImgPtr(new DiskImageView<PixelT>(fileName) ) );
	}catch(std::exception const& e){
	  freeup_handles_not_thread_safe();
	  m_fileHandles.push_back( ImgPtr() );
	}
      }else{
        m_fileHandles.push_back( ImgPtr() );
      }
      
      int index = m_fileHandles.size() - 1;
      m_used[index] = 0; // open, but not in active use yet
    }

    // Get a handle to the current image. If too many handles are
    // open, before opening a new one close the ones whose pixel boxes
    // are furthest than the current desired box we want to use the
    // handle.  This way, we ensure we don't close handles more often
    // than what we have to.
    DiskImageView<PixelT> & get_handle(int index, BBox2 const& bbox){

      Mutex::WriteLock write_lock(m_mutex);

      m_used[index]++;
      
      // Just return the handle if it is present
      if (m_fileHandles[index].get() != NULL)
	return *m_fileHandles[index].get();

      double min_pos_dist = 1.0;
      
      // Find the image whose box is furthest form the current one.
      std::vector<double> boxes_dist;
      for (size_t b = 0; b < m_imgBoxes.size(); b++) {
        
	// If this image is being used by another thread, can't close it.
	if (m_used[b] > 0) {
	  boxes_dist.push_back(0);
	  continue;
	}

	// If boxes intersect, they are of course very close.
	// We'd like to pick such an image only if we really
	// can't pick other ones whose boxes don't intersect,
	// hence 1.5*min_pos_dist below.
	BBox2 curr_box     = m_imgBoxes[b];
	BBox2 intersection = curr_box;
	intersection.crop(bbox);
	if (!intersection.empty()) {
	  boxes_dist.push_back(1.5*min_pos_dist);
	  continue;
	}

	// If the current handle is not open, can't close it.
	// Notice we give it a positive distance, but a small
	// one, because we want handles that are closed to be
	// picked only if there are no open handles. 
	// That is why we add 2*min_pos_dist below.
	if (m_fileHandles[b].get() == NULL){
	  boxes_dist.push_back(min_pos_dist);
	  continue;
	}
	double dist
	  = norm_2((curr_box.min() + curr_box.max())/2.0
		   - (bbox.min() + bbox.max())/2.0)
	  + 2*min_pos_dist;
	if (curr_box.empty() || bbox.empty()) 
	  dist = 2*min_pos_dist; // overabundance of caution
	
	boxes_dist.push_back(dist);
      }

      if (boxes_dist.size() != m_fileHandles.size() ) 
	vw_throw(ArgumentErr() << "Book-keeping failure in " <<
		 __FILE__ << "\n");
      
      std::vector<double>::iterator max_it
	= std::max_element(boxes_dist.begin(), boxes_dist.end());
      int max_index = std::distance(boxes_dist.begin(), max_it);

      double max_val = boxes_dist[max_index];
      if (max_val == 0) 
	vw_throw(ArgumentErr() << "Could not find a handle to "
		 << "free up. Too many input images?\n");
      
      // This one better be open
      if (m_fileHandles[max_index].get() == NULL)
	vw_throw(ArgumentErr() << "Expecting an open handle.\n");

      // Close it
      vw_out() << "\nClosing: " << m_fileNames[max_index] << std::endl;
      m_fileHandles[max_index] = ImgPtr();

      // Open the current handle
      vw_out() << "\nLoading: " << m_fileNames[index] << std::endl;
      m_fileHandles[index] = ImgPtr(new DiskImageView<PixelT>(m_fileNames[index]));

      return *m_fileHandles[index].get();
    }

    std::string get_file_name(int index){
      if (index < 0 || index >= (int)m_fileNames.size()) {
	vw_throw(ArgumentErr() << "Could not find file with index " << index << "\n");
      }
      return m_fileNames[index];
    }
    
    // Declare a given handle as unused
    void release(int index){
      Mutex::WriteLock write_lock(m_mutex);
      std::map<int, int>::iterator it = m_used.find(index);
      if (it == m_used.end()) return; // Not in the set

      if (m_used[index] > 0) 
	m_used[index]--; // one less used
    }

  }; // End class DiskImageManager
  

} // namespace vw

#endif  // __VW_FILEIO_DISKIMAGEMANAGER_H__
