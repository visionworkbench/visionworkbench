#ifdef _MSC_VER
#pragma warning(disable:4244)
#pragma warning(disable:4267)
#pragma warning(disable:4996)
#endif

#include <iostream>
#include <vector>

#include <boost/smart_ptr.hpp>

#include <vw/Image/Algorithms.h>
#include <vw/FileIO/DiskImageResource.h>
#include <vw/Mosaic/ImageComposite.h>

void vw::mosaic::ImageComposite::generate_masks() const {
  vw_out(InfoMessage) << "Generating masks..." << std::endl;
  std::vector<boost::shared_ptr<ImageView<float> > > grassfire_im( sources.size() );
  for( unsigned p1=0; p1<sources.size(); ++p1 ) {
    if( ! grassfire_im[p1] ) {
      ImageView<PixelRGBA<float> > image = *sources[p1];
      grassfire_im[p1].reset( new ImageView<float>( grassfire( select_channel( image, 3 ) ) ) );
    }
    ImageView<float> mask = copy( *grassfire_im[p1] );
    for( unsigned p2=0; p2<sources.size(); ++p2 ) {
      if( p1 == p2 ) continue;
      int ox = bboxes[p2].min().x() - bboxes[p1].min().x();
      int oy = bboxes[p2].max().y() - bboxes[p1].max().y();
      if( ox >= bboxes[p1].width() ||
          oy >= bboxes[p1].height() ||
          -ox >= bboxes[p2].width() ||
          -oy >= bboxes[p2].height() ) {
        grassfire_im[p2].reset();
      }
      else {
        if( ! grassfire_im[p2] ) {
          ImageView<PixelRGBA<float> > image = *sources[p2];
          grassfire_im[p2].reset( new ImageView<float>( grassfire( select_channel( image, 3 ) ) ) );
        }
        int left = std::max( ox, 0 );
        int top = std::max( oy, 0 );
        int right = std::min( bboxes[p2].width()+ox, bboxes[p1].width() );
        int bottom = std::min( bboxes[p2].height()+oy, bboxes[p1].height() );
        for( int j=top; j<bottom; ++j ) {
          for( int i=left; i<right; ++i ) {
            if( ( (*grassfire_im[p2])(i-ox,j-oy) > mask(i,j) ) ||
                ( (*grassfire_im[p2])(i-ox,j-oy) == mask(i,j) && p2 > p1 ) )
              mask(i,j) = 0;
          }
        }
      }
    }
    mask = threshold( mask );
    std::ostringstream filename;
    filename << "mask." << p1 << ".png";
    write_image( filename.str(), mask );
  }
}


boost::shared_ptr<vw::mosaic::ImageComposite::Pyramid> vw::mosaic::ImageComposite::PyramidGenerator::generate() const {
  vw_out(DebugMessage) << "ImageComposite generating pyramid " << m_index << std::endl;
  boost::shared_ptr<Pyramid> ptr( new Pyramid );
  ImageView<PixelRGBA<float> > source = copy(*m_composite.sources[m_index]);
  m_composite.sources[m_index].deprioritize();
  PositionedImage<PixelRGBA<float> > image_high( m_composite.bbox.width(), m_composite.bbox.height(), source, m_composite.bboxes[m_index] );
  PositionedImage<PixelRGBA<float> > image_low = image_high.reduce();
  ImageView<float> mask_image;
  
  std::ostringstream mask_filename;
  mask_filename << "mask." << m_index << ".png";
  read_image( mask_image, mask_filename.str() );
  PositionedImage<float> mask( m_composite.bbox.width(), m_composite.bbox.height(), mask_image, m_composite.bboxes[m_index] );
  
  for( int l=0; l<m_composite.levels; ++l ) {
    PositionedImage<PixelRGBA<float> > diff = image_high;
    if( l > 0 ) mask = mask.reduce();
    if( l < m_composite.levels-1 ) {
      PositionedImage<PixelRGBA<float> > next_image_low = image_low.reduce();
      image_low.unpremultiply();
      diff.subtract_expanded( image_low );
      image_high = image_low;
      image_low = next_image_low;
    }
    diff *= mask;
    ptr->images.push_back( diff );
    ptr->masks.push_back( mask );
  }
  return ptr;
}


void vw::mosaic::ImageComposite::insert( ImageViewRef<PixelRGBA<float> > const& image, int x, int y ) {
  vw_out(VerboseDebugMessage) << "ImageComposite inserting image " << pyramids.size() << std::endl;
  sources.push_back( m_cache.insert( SourceGenerator( image ) ) );
  alphas.push_back( m_cache.insert( AlphaGenerator( *this, pyramids.size() ) ) );
  pyramids.push_back( m_cache.insert( PyramidGenerator( *this, pyramids.size() ) ) );
  int cols = image.cols(), rows = image.rows();
  BBox<int,2> image_bbox( Vector<int,2>(x, y), Vector<int,2>(x+cols, y+rows) );
  bboxes.push_back( image_bbox );
  if( bboxes.size() == 1 ) {
    bbox = bboxes.back();
    mindim = std::min(cols,rows);
  }
  else {
    bbox.grow( image_bbox );
    mindim = std::min( mindim, std::min(cols,rows) );
  }
}


void vw::mosaic::ImageComposite::prepare() {
  // Translate bboxes to origin
  for( unsigned i=0; i<sources.size(); ++i )
    bboxes[i] -= bbox.min();

  levels = (int) floorf( log( mindim/2.0 ) / log(2.0) ) - 1;

  if( !m_draft_mode ) {
    generate_masks();
  }
}


// Suppose a destination image patch at a given level of the pyramid
// has a bounding box that begins at offset x and has width w.  It
// is affected by a range of pixels at the next level of the pyramid
// starting at x/2 with width (x+w)/2-x/2+1 = (w+x%2)/2+1.  This in
// turn is affected by source image pixels at the current level in
// the range starting at 2*(x/2)-1 = x-x%2-1 with width
// (2*(x+w)/2+1)-(2*(x/2)-1)+1 = w-(x+w)%2+x%2+3.

// Generates a full-resolution patch of the mosaic corresponding
// to the given bounding box.
vw::ImageView<vw::PixelRGBA<float> > vw::mosaic::ImageComposite::blend_patch( BBox<int,2> const& patch_bbox ) const {
  vw_out(DebugMessage) << "ImageComposite compositing patch " << patch_bbox << "..." << std::endl;
  // Compute bboxes and allocate the pyramids
  std::vector<BBox<int,2> > bbox_pyr;
  std::vector<ImageView<PixelRGBA<float> > > sum_pyr(levels);
  std::vector<ImageView<float> > msum_pyr(levels);
  for( int l=0; l<levels; ++l ) {
    if( l==0 ) bbox_pyr.push_back( patch_bbox );
    else bbox_pyr.push_back( BBox<int,2>( Vector<int,2>( bbox_pyr[l-1].min().x() / 2,
                                                         bbox_pyr[l-1].min().y() / 2 ),
                                          Vector<int,2>( bbox_pyr[l-1].min().x() / 2 + ( bbox_pyr[l-1].width() + bbox_pyr[l-1].min().x() % 2 ) / 2 + 1,
                                                         bbox_pyr[l-1].min().y() / 2 + ( bbox_pyr[l-1].height() + bbox_pyr[l-1].min().y() % 2 ) / 2 + 1) ) );
    sum_pyr[l] = ImageView<PixelRGBA<float> >( bbox_pyr[l].width(), bbox_pyr[l].height() );
    msum_pyr[l] = ImageView<float>( bbox_pyr[l].width(), bbox_pyr[l].height() );
  }
  
  // Compute the bounding box for source pixels that could 
  // impact the patch.
  BBox<int,2> padded_bbox = patch_bbox;
  for( int l=0; l<levels-1; ++l ) {
    padded_bbox.min().x() = padded_bbox.min().x()/2;
    padded_bbox.min().y() = padded_bbox.min().y()/2;
    padded_bbox.max().x() = padded_bbox.max().x()/2+1;
    padded_bbox.max().y() = padded_bbox.max().y()/2+1;
  }
  for( int l=0; l<levels-1; ++l ) {
    padded_bbox.min().x() = 2*padded_bbox.min().x()-1;
    padded_bbox.min().y() = 2*padded_bbox.min().y()-1;
    padded_bbox.max().x() = 2*padded_bbox.max().x();
    padded_bbox.max().y() = 2*padded_bbox.max().y();
  }
  
  // Make a list of the images whose bounding boxes permit them to
  // impact the patch, prioritizing ones that are already in memory.
  std::list<unsigned> image_list;
  for( unsigned p=0; p<sources.size(); ++p ) {
    if( ! padded_bbox.intersects( bboxes[p] ) ) continue;
    if( ! pyramids[p].valid() ) image_list.push_back( p );
    else image_list.push_front( p );
  }

  // Add each source image pyramid to the blend pyramid.
  std::list<unsigned>::iterator ili=image_list.begin(), ilend=image_list.end();
  for( ; ili!=ilend; ++ili ) {
    unsigned p = *ili;
    boost::shared_ptr<Pyramid> pyr = pyramids[p];
    for( int l=0; l<levels; ++l ) {
      pyr->images[l].addto( sum_pyr[l], bbox_pyr[l].min().x(), bbox_pyr[l].min().y() );
      pyr->masks[l].addto( msum_pyr[l], bbox_pyr[l].min().x(), bbox_pyr[l].min().y() );
    }
  }

  // Collapse the pyramid
  ImageView<PixelRGBA<float> > composite( sum_pyr[levels-1].cols(), sum_pyr[levels-1].rows() );
  for( int l=levels; l; --l ) {
    if( l < levels ) {
      composite = ImageView<PixelRGBA<float> >( crop( resample( composite, 2 ), 
                                                      bbox_pyr[l-1].min().x()-2*bbox_pyr[l].min().x(), 
                                                      bbox_pyr[l-1].min().y()-2*bbox_pyr[l].min().y(), 
                                                      sum_pyr[l-1].cols(), sum_pyr[l-1].rows() ) );
    }
    composite += sum_pyr[l-1] / msum_pyr[l-1];
    sum_pyr.pop_back();
    msum_pyr.pop_back();
  }

  // Trim to the maximal source alpha, reloading images if needed
  ImageView<float> alpha( patch_bbox.width(), patch_bbox.height() );
  for( unsigned p=0; p<sources.size(); ++p ) {
    if( ! patch_bbox.intersects( bboxes[p] ) ) continue;
    
    ImageView<float> source_alpha = *alphas[p];
    
    BBox<int,2> overlap = patch_bbox;
    overlap.crop( bboxes[p] );
    for( int j=0; j<overlap.height(); ++j ) {
      for( int i=0; i<overlap.width(); ++i ) {
        if( source_alpha( overlap.min().x()+i-bboxes[p].min().x(), overlap.min().y()+j-bboxes[p].min().y() ) > 
            alpha( overlap.min().x()+i-patch_bbox.min().x(), overlap.min().y()+j-patch_bbox.min().y() ) ) {
          alpha( overlap.min().x()+i-patch_bbox.min().x(), overlap.min().y()+j-patch_bbox.min().y() ) =
            source_alpha( overlap.min().x()+i-bboxes[p].min().x(), overlap.min().y()+j-bboxes[p].min().y() );
        }
      }
    }
  }
  composite *= alpha / select_channel( composite, 3 );

  return composite;
}


// Generates a full-resolution patch of the mosaic corresponding
// to the given bounding box WITHOUT blending.
vw::ImageView<vw::PixelRGBA<float> > vw::mosaic::ImageComposite::draft_patch( BBox<int,2> const& patch_bbox ) const {
  vw_out(DebugMessage) << "ImageComposite compositing patch " << patch_bbox << "..." << std::endl;
  ImageView<PixelRGBA<float> > composite(patch_bbox.width(),patch_bbox.height());

  // Add each image to the composite.
  for( unsigned p=0; p<sources.size(); ++p ) {
    if( ! patch_bbox.intersects( bboxes[p] ) ) continue;
    PositionedImage<PixelRGBA<float> > image( bbox.width(), bbox.height(), *sources[p], bboxes[p] );
    image.addto( composite, patch_bbox.min().x(), patch_bbox.min().y(), true );
  }

  return composite;
}
