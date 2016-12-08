// __BEGIN_LICENSE__
//  Copyright (c) 2009-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__


/// \file BlobIndex.h
///

#ifndef __VW_IMAGE_BLOB_INDEX_THREADED_H__
#define __VW_IMAGE_BLOB_INDEX_THREADED_H__

// VW
#include <vw/Core/Log.h>
#include <vw/Core/Thread.h>
#include <vw/Core/ThreadPool.h>
#include <vw/Core/Stopwatch.h>
#include <vw/Image/AlgorithmFunctions.h>
#include <vw/Image/Manipulation.h>
#include <vw/Image/ImageView.h>
#include <vw/Image/PixelAccessors.h>
#include <vw/Image/PixelMask.h>

// Standard
#include <vector>
#include <deque>
#include <list>
#include <ostream>

// Boost
#include <boost/noncopyable.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_selectors.hpp>

// BlobIndex (Multi) Threaded
///////////////////////////////////////

// This code is a custom write up of the blob_index function from VW. It has several
// different features.
// --> Multithread (as expected)
// --> Lower Memory Impact
//     via a new internal compressed format
// --> Allows for limiting on size.

namespace vw {

namespace blob {

  // Blob Compressed
  ////////////////////////////////////
  // A nice way to describe a blob,
  // but reducing our memory foot print
  class BlobCompressed {
    // This describes a blob as lines of rows
    // to reduce the memory foot print
    Vector2i m_min;
    std::vector< std::list<int32> > m_row_start; // assumed to be ordered
    std::vector< std::list<int32> > m_row_end;

    void shift_x ( int32 const& value );
    void refactor();

  public:
    BlobCompressed( Vector2i const& top_left,
                    std::vector<std::list<int32> > const& row_start,
                    std::vector<std::list<int32> > const& row_end );
    BlobCompressed();

    // Standard Access point
    Vector2i const& min() const;
    Vector2i      & min();

    std::list<int32> const& start( uint32 const& index ) const;
    std::list<int32> const& end  ( uint32 const& index ) const;

    int32 num_rows() const;
    int32 size    () const; // Please use sparingly

    BBox2i bounding_box() const;
    bool intersects( BBox2i const& input ) const;

    // Specific conditionals used by BlobIndexThreaded
    bool is_on_right ( BlobCompressed const& right  ) const;
    bool is_on_bottom( BlobCompressed const& bottom ) const;

    // Append a row (since these guys are built a row at a time )
    void add_row( Vector2i const& start, int const& width );
    // Use to expand this blob into a non overlapped area
    void absorb( BlobCompressed const& victim );
    // Dump listing of every pixel used
    void decompress( std::list<Vector2i>& output ) const;
    // Print internal data
    void print() const;
  };

  // Blob Index Custom
  ////////////////////////////////////
  /// A different version of Blob index that uses the compressed format and the new options
  class BlobIndex {
    std::vector<BlobCompressed> m_c_blob;
    uint m_blob_count;

  public:
    // Constructor performs processing
    template <class SourceT>
    BlobIndex( ImageViewBase<SourceT> const& src, // Masked input image
               ImageView<uint32>           & dst, // Output image contains labeled blobs
               uint /*max_area*/=0 ) {

      if ( src.impl().planes() > 1 )
        vw_throw( NoImplErr()
                  << "Blob index currently only works with 2D images." );
      dst.set_size( src.impl().cols(),
                    src.impl().rows() );
      fill(dst,0);

      // Initialize Graph used to pair blobs
      typedef boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS> Graph;
      Graph connections;
      m_blob_count=1;

      { // Initial Pass
        // Leading and trailing iterators for the input and output images.
        typename SourceT::pixel_accessor   s_acc = src.impl().origin();
        typename SourceT::pixel_accessor p_s_acc = src.impl().origin(); // previous
        typename ImageView<uint32>::pixel_accessor   d_acc = dst.origin();
        typename ImageView<uint32>::pixel_accessor p_d_acc = dst.origin(); // previous

        // Top Corner
        if ( is_valid(*s_acc) ) {
          *d_acc = m_blob_count;
          m_blob_count++;
        }

        // Top Row
        s_acc.next_col();
        d_acc.next_col();
        for ( int32 i = 1; i < dst.cols(); i++ ) {
          if ( is_valid(*s_acc) ) {
            if ( is_valid(*p_s_acc) ) { // Last pixel and this pixel valid
              *d_acc = *p_d_acc; // Extend the current blob label.
            } else { // Start of a new blob
              *d_acc = m_blob_count;
              m_blob_count++;
            }
          }
          s_acc.next_col(); // Advance all iterators
          d_acc.next_col();
          p_s_acc.next_col();
          p_d_acc.next_col();
        }
      }

      { // Everything else (9 connected)
        typename SourceT::pixel_accessor           s_acc_row = src.impl().origin();
        typename ImageView<uint32>::pixel_accessor d_acc_row = dst.origin();
        s_acc_row.advance(0,1); // We already did the first row, so go to the second.
        d_acc_row.advance(0,1);

        // Loop through the rows
        for (int j = dst.rows()-1; j; --j ) {
          // Leading and trailing iterators for the input and output images.
          typename SourceT::pixel_accessor             s_acc = s_acc_row;
          typename SourceT::pixel_accessor           p_s_acc = s_acc_row;
          typename ImageView<uint32>::pixel_accessor   d_acc = d_acc_row;
          typename ImageView<uint32>::pixel_accessor p_d_acc = d_acc_row;

          // Process
          
          // Loop through columns
          for ( int i = dst.cols(); i; --i ) {
          
            if ( is_valid(*s_acc) ) { // Current pixel is valid
            
              // TODO: Can we make this more efficient?
            
              if ( i != dst.cols() ) { // Not the first column
                // Check the pixel to the left
                p_s_acc.advance(-1,0);
                p_d_acc.advance(-1,0);
                if ( is_valid(*p_s_acc) ) { // It is valid...
                  if ( (*d_acc != 0) && (*d_acc != *p_d_acc) ) // Will this case ever trigger? How is d_acc set here?
                    boost::add_edge(*p_d_acc,*d_acc,connections);
                  else
                    *d_acc = *p_d_acc;
                }
                // Check the pixel to the upper left
                p_s_acc.advance(0,-1);
                p_d_acc.advance(0,-1);
                if ( is_valid(*p_s_acc) ) { // If valid...
                  if ( (*d_acc != 0) && (*d_acc != *p_d_acc) ) // If current pixel and prev pixel have different labels
                    boost::add_edge(*p_d_acc,*d_acc,connections); // Add a connection between the labels in boost
                  else
                    *d_acc = *p_d_acc; // Copy the label from the previous pixel
                }
              } else { // This is the first column
                p_s_acc.advance(-1,-1); // Move to the non-existant upper left location so that
                p_d_acc.advance(-1,-1); //  the following commands to to the correct locations.
              }
              
              // Check the upper pixel
              p_s_acc.advance(1,0);
              p_d_acc.advance(1,0);
              if ( is_valid(*p_s_acc) ) {
                if ( (*d_acc != 0) && (*d_acc != *p_d_acc) )
                  boost::add_edge(*p_d_acc,*d_acc,connections);
                else
                  *d_acc = *p_d_acc;
              }
              
              // Check the upper right pixel
              p_s_acc.advance(1,0);
              p_d_acc.advance(1,0);
              if ( i != 1 ) { // If not on the last column
                if ( is_valid(*p_s_acc) ) {
                  if ( (*d_acc != 0) && (*d_acc != *p_d_acc) )
                    boost::add_edge(*p_d_acc,*d_acc,connections);
                  else
                    *d_acc = *p_d_acc;
                }
              }
              
              // Move the trailing iterators up to the current pixel location.
              p_s_acc.advance(-1,1);
              p_d_acc.advance(-1,1);
              if ( *d_acc == 0 ) { // If we have not yet assigned this pixel, start a new blob.
                *d_acc = m_blob_count;
                m_blob_count++;
              }
            } // End case where current pixel is valid
            
            // Advance to the next column
            s_acc.next_col();
            p_s_acc.next_col();
            d_acc.next_col();
            p_d_acc.next_col();
          } // end column loop

          // Anvance to the next row
          s_acc_row.next_row();
          d_acc_row.next_row();
        } // End row loop
      } // End phantom bracket

      // Making sure connections has vertices for all indexes made
      add_edge(m_blob_count-1,m_blob_count-1,connections);
      std::vector<uint32> component(boost::num_vertices(connections));
      m_blob_count = boost::connected_components(connections, &component[0])-1;
      m_c_blob.resize(m_blob_count);
      // component contains the true (consolidated) labels for each pixel
      // m_c_blob will store the output blobs.

      { // Update index map to optimal numbering
        ImageView<uint32>::pixel_accessor p_d_acc = dst.origin(); // previous

        // Loop through rows
        for ( int32 r = 0; r < dst.rows(); r++ ) {
          ImageView<uint32>::pixel_accessor d_acc = p_d_acc;
          bool building_segment=false;
           int32 start_col = 0;
          uint32 index     = 0;
          
          // Loop through columns
          // - In each column, build up a run-length-encoding for each blob.
          for ( int32 c = 0; c < dst.cols(); c++ ) {
            if ( (*d_acc) != 0 ) { // If labeled pixel...
              if ( building_segment && (index != component[*d_acc]) ) {
                vw_throw(LogicErr() << "Insert seems wrong.\n");
                // I believe the way it's processed this shouldn't happen
              } else if (!building_segment) {
                // Start recording a new RLE segment
                building_segment = true;
                index     = component[*d_acc]; // Get the true label for this pixel
                start_col = c;
              }
            } else if ( building_segment ) {
              // Finished with this blob for now, add an RLE segment
              building_segment = false;
              m_c_blob[index-1].add_row( Vector2i(start_col,r), c-start_col );
            }
            d_acc.next_col();
          }
          
          if ( building_segment ) { // Hit the end of a row building a blob, cut it off.
            m_c_blob[index-1].add_row( Vector2i(start_col,r),
                                       dst.cols()-start_col );
          }
          p_d_acc.next_row();
        }
      }

    }

    // Access points to intersting information
    uint32 num_blobs() const;

    // Access to blobs
    BlobCompressed const& blob( uint32 const& index ) const;
  };



  // Blob Index Task
  /////////////////////////////////////
  /// A task wrapper to allow threading
  template <class SourceT>
  class BlobIndexTask : public Task, private boost::noncopyable {

    ImageViewBase<SourceT> const& m_view;
    BBox2i const& m_bbox;
    Mutex&        m_append_mutex;

    std::deque<BlobCompressed> &m_c_blob; // reference to global
    std::deque<BBox2i>     &m_blob_bbox;
    int m_id;
    int m_max_area;
  public:
    BlobIndexTask( ImageViewBase<SourceT> const& view,
                   BBox2i const& bbox, Mutex &mutex,
                   std::deque<BlobCompressed> & blobs,
                   std::deque<BBox2i> & blob_boxes,
                   int const& id, int const& max_area ) :
    m_view(view), m_bbox(bbox), m_append_mutex(mutex),
      m_c_blob(blobs), m_blob_bbox(blob_boxes), m_id(id), m_max_area(max_area) {}

    void operator()() {
      Stopwatch sw;
      sw.start();
      ImageView<uint32> index_image(m_bbox.width(),
                                            m_bbox.height() );

      // Render so threads don't wait on each other
      ImageView<typename SourceT::pixel_type> cropped_copy = crop(m_view,m_bbox);
      // Decided only to do trimming in the global perspective. This
      // avoids weird edge effects.
      BlobIndex bindex( cropped_copy, index_image);

      // Build local bboxes
      std::vector<BBox2i> local_bboxes( bindex.num_blobs() );
      for ( uint32 i = 0; i < bindex.num_blobs(); i++ ) {
        local_bboxes[i] = bindex.blob(i).bounding_box() + m_bbox.min();
      }

      { // Append results (single thread)
        Mutex::Lock lock(m_append_mutex);
        for ( uint i = 0; i < local_bboxes.size(); i++ ) {
          m_c_blob.push_back( bindex.blob(i) );
          m_c_blob.back().min() += m_bbox.min(); // Fix offset
          m_blob_bbox.push_back( local_bboxes[i] );
        }
      }

      sw.stop();
      vw_out(VerboseDebugMessage,"inpaint") << "Task " << m_id << ": finished, " << sw.elapsed_seconds() << "s\n";
    }
  };
} // end namespace blob


  // Simple interface
  template <class SourceT>
  ImageView<uint32> blob_index( ImageViewBase<SourceT> const& src ) {
    ImageView<uint32> result( src.impl().cols(), src.impl().rows() );
    blob::BlobIndex( src, result );
    return result;
  }

// Blob Index Threaded
///////////////////////////////////
/// Performs Blob Index using all threads and a minimal amount of memory
class BlobIndexThreaded {

  std::deque<BBox2i>           m_blob_bbox;
  std::deque<blob::BlobCompressed> m_c_blob;

  Mutex m_insert_mutex;
  int m_max_area;
  int m_tile_size;

  // Tasks might section a blob in half.
  // This will match them
  void consolidate( Vector2i const& image_size,
                    Vector2i const& proc_block_size );

 public:
 
  /// Constructor does most of the processing work
  /// - This is the function to call to detect blobs!
  /// - Blobs larger than max_area (if > zero) are discarded.
  template <class SourceT>
  BlobIndexThreaded( ImageViewBase<SourceT> const& src,
                     int32 const& max_area    = 0,
                     int32 const& tile_size   = vw_settings().default_tile_size(),
                     int32 const& num_threads = vw_settings().default_num_threads()
                     )
    : m_max_area(max_area), m_tile_size(tile_size) {
    
    std::vector<BBox2i> bboxes = subdivide_bbox( src.impl(), m_tile_size, m_tile_size );
    // User needs to remember to give a pixel mask'd input
    typedef blob::BlobIndexTask<SourceT> task_type;
    if (bboxes.size() > 1){
      Stopwatch sw;
      sw.start();
      FifoWorkQueue queue(num_threads);
      
      for ( size_t i = 0; i < bboxes.size(); ++i ) {
        boost::shared_ptr<task_type> task(new task_type(src, bboxes[i],
                                                        m_insert_mutex,
                                                        m_c_blob, m_blob_bbox,
                                                        i, m_max_area ));
        queue.add_task(task);
      }
      queue.join_all();
      
      sw.stop();
      vw_out(DebugMessage,"inpaint") << "Blob detection took " << sw.elapsed_seconds() << "s\n";
      consolidate( Vector2i( src.impl().cols(), src.impl().rows() ),
                   Vector2i( m_tile_size, m_tile_size ) );
    }else if (bboxes.size() == 1){
      // This is a special case when we want to fill in holes
      // in a single small image.
      boost::shared_ptr<task_type> task(new task_type(src, bboxes[0],
                                                      m_insert_mutex,
                                                      m_c_blob, m_blob_bbox,
                                                      0, m_max_area ));
      (*task.get())();
    }
    

    // Cull blobs that are too big.
    if ( m_max_area > 0 ) {
      for ( std::deque<blob::BlobCompressed>::iterator iter = m_c_blob.begin();
            iter != m_c_blob.end(); iter++ ) {
        if ( iter->size() > m_max_area ) {
          iter = m_c_blob.erase( iter );
          iter--;
        }
      }
    }
  }

  /// Wipe blobs bigger than this size.
  void wipe_big_blobs(int max_size){

    // Keep the bounding boxes up-to-date.
    m_blob_bbox.clear();
    
    for ( std::deque<blob::BlobCompressed>::iterator iter = m_c_blob.begin();
          iter != m_c_blob.end(); iter++ ) {

      BBox2i blob_bbox = iter->bounding_box();
      
      if ( blob_bbox.width() > max_size ||
           blob_bbox.height() > max_size) {
        iter = m_c_blob.erase( iter );
        iter--;
      }else{
        m_blob_bbox.push_back(blob_bbox);
      }
      
    }
    
  }
  
  // Access for the users
  uint32 num_blobs() const;
  /// ?
  void blob( uint32 const& index,
             std::list<Vector2i>& output ) const;
  /// ?
  blob::BlobCompressed const& compressed_blob( uint32 const& index ) const;

  typedef std::deque<blob::BlobCompressed>::iterator             blob_iterator;
  typedef std::deque<blob::BlobCompressed>::const_iterator const_blob_iterator;
        blob_iterator begin();
  const_blob_iterator begin() const;
        blob_iterator end();
  const_blob_iterator end() const;

  BBox2i const& blob_bbox( uint32 const& index ) const;
  typedef std::deque<BBox2i>::iterator             bbox_iterator;
  typedef std::deque<BBox2i>::const_iterator const_bbox_iterator;
        bbox_iterator bbox_begin();
  const_bbox_iterator bbox_begin() const;
        bbox_iterator bbox_end();
  const_bbox_iterator bbox_end() const;
};



/// From a masked image, generate an image where each pixel has a value
/// equal to the size of the blob that contains it (up to a size limit).
template <class ImageT>
class BlobSizesView: public ImageViewBase<BlobSizesView<ImageT> >{
  ImageT const& m_input_image;
  int    m_expand_size; ///< Tile expansion used to more accurately size blobs.
  uint32 m_size_limit;  ///< Cap the size value written in output pixels.
public:
  /// Constructor
  /// - "expand_size" is the distance searched out from each tile.  A larger value
  ///   will increase the run time but will improve accuracy long, narrow blobs
  ///   near tile borders.
  /// - "size_limit" is the maximum size score that can be assigned to a blob.
  BlobSizesView( ImageViewBase<ImageT> const& img, int expand_size, 
                 uint32 size_limit = std::numeric_limits<uint32>::max()):
    m_input_image(img.impl()), m_expand_size(expand_size), m_size_limit(size_limit){}

  // Image View interface
  typedef uint32 pixel_type;
  typedef pixel_type      result_type;
  typedef ProceduralPixelAccessor<BlobSizesView> pixel_accessor;

  inline int32 cols  () const { return m_input_image.cols(); }
  inline int32 rows  () const { return m_input_image.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  inline pixel_type operator()( double i, double j, int32 p = 0 ) const {
    vw_throw(NoImplErr() << "operator()(...) is not implemented");
    return pixel_type();
  }

  typedef CropView<ImageView<pixel_type> > prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i const& bbox) const {

    ImageView<pixel_type> output_tile(bbox.width(), bbox.height());

    // Rasterize the region with an expanded tile size so that we can count
    //  blob sizes that extend some distance outside the tile borders.
    BBox2i big_bbox = bbox;
    big_bbox.expand(m_expand_size);
    big_bbox.crop(bounding_box(m_input_image));
    ImageView<typename ImageT::pixel_type> input_tile = crop(m_input_image, big_bbox);
    //std::cout << "Searching for blobs in " << big_bbox << std::endl;

    // Use a single thread to search for blobs in this image with no max blob size.
    int tile_size = std::max(big_bbox.width(), big_bbox.height());
    BlobIndexThreaded blob_index(input_tile, 0, tile_size, 1);
    
    // Loop through all the blobs and assign pixel scores
    BlobIndexThreaded::const_blob_iterator blob_iter = blob_index.begin();
    //std::cout << "Found " << blob_index.num_blobs() << " blobs.\n";
    while (blob_iter != blob_index.end()) { // Loop through blobs
      uint32 blob_size = blob_iter->size();
      //std::cout << "Blob size =  " << blob_size << "\n";
      if (blob_size > m_size_limit)
        blob_size = m_size_limit;
        
      // Loop through rows in blobs
      int num_rows = blob_iter->num_rows();
      //std::cout << "num_rows = " << num_rows << std::endl;
      std::list<int32>::const_iterator start_iter, stop_iter;
      for (int r=0; r<num_rows; ++r) {
      
        // Loop through sections in row
        int row = r + blob_iter->min()[1] + big_bbox.min()[1]; // Absolute row
        start_iter = blob_iter->start(r).begin();
        stop_iter  = blob_iter->end(r).begin();
        while (start_iter != blob_iter->start(r).end()) { 
          //std::cout << "start = " << *start_iter << ", stop = " << *stop_iter << std::endl;
          
          // Loop through pixels in section
          for (int c=*start_iter; c<*stop_iter; ++c) { 
            int col = c + blob_iter->min()[0] + big_bbox.min()[0]; // Absolute col
            Vector2i pixel(col,row);
            //std::cout << "Searching for blobs in " << big_bbox << std::endl;
            //std::cout << "Pixel = " << pixel << ", raw = " << Vector2i(c,r) << std::endl;
            if (!bbox.contains(pixel))
              continue; // Skip blob pixels outside the current tile
            Vector2i tile_pixel = pixel - bbox.min();
            //std::cout << "Pixel = " << pixel << ", tile_pixel = " << tile_pixel << ", raw = " << Vector2i(c,r) << std::endl;
            output_tile(tile_pixel[0], tile_pixel[1]) = blob_size; // Set the size of the pixel!
          } // End loop through pixels
          
          ++start_iter;
          ++stop_iter;
        } // End loop through sections
        
      } // End loop through rows
      
      ++blob_iter;
    } // End loop through blobs

    // Perform tile size faking trick to make small tile look the size of the entire image
    return prerasterize_type(output_tile,
                             -bbox.min().x(), -bbox.min().y(),
                             cols(), rows() );
  }

  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
}; // End class BlobSizesView

template <class ImageT>
BlobSizesView<ImageT>
get_blob_sizes(ImageViewBase<ImageT> const& image, int expand_size, uint32 size_limit) {
  return BlobSizesView<ImageT>(image.impl(), expand_size, size_limit);
}




} // end namespace vw

#endif//__BLOB_INDEX_THREADED_H__
