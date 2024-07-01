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

// Turn off warnings about things we can't control
#include <vw/Math/FLANNTree.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <flann/flann.hpp>
#pragma GCC diagnostic pop

namespace vw {
namespace math {

// This variable and lock will be used to print a message once
bool flann_message_printed = false;
Mutex flann_mutex;

//=============================================================
// - Only the float and double versions of this code are listed here.
// - This multiplies the code size, but hides the messy flann headers from the .h file.

// Helper functions to clean up the pointer casting below
flann::Index<flann::L2<float>>* cast_index_ptr_L2_f(void* void_ptr) {
  return reinterpret_cast<flann::Index<flann::L2<float>>*>(void_ptr);
}
const flann::Index<flann::L2<float>>* cast_index_ptr_L2_f(const void* void_ptr) {
  return reinterpret_cast<const flann::Index<flann::L2<float>>*>(void_ptr);
}

// Helper functions to clean up the pointer casting below
flann::Index<flann::L2<double>>* cast_index_ptr_L2_d(void* void_ptr) {
  return reinterpret_cast<flann::Index<flann::L2<double>>*>(void_ptr);
}
const flann::Index<flann::L2<double>>* cast_index_ptr_L2_d(const void* void_ptr) {
  return reinterpret_cast<const flann::Index<flann::L2<double>>*>(void_ptr);
}

// Helper functions to clean up the pointer casting below
flann::Index<flann::Hamming<unsigned char>>* cast_index_ptr_HAMM_u(void* void_ptr) {
  return reinterpret_cast<flann::Index<flann::Hamming<unsigned char>>*>(void_ptr);
}

const flann::Index<flann::Hamming<unsigned char>>* cast_index_ptr_HAMM_u(const void* void_ptr) {
  return reinterpret_cast<const flann::Index<flann::Hamming<unsigned char>>*>(void_ptr);
}

template <>
size_t FLANNTree<float>::knn_search_help(void* data_ptr, size_t rows, size_t cols,
                                         Vector<int>& indices,
                                         Vector<double>& dists,
                                         size_t knn) {
  // Constrain the number of results that we can return to the number of loaded objects
  size_t maxNumReturns = m_features_cast.rows();
  if (knn > maxNumReturns)
    knn = maxNumReturns;

  if (indices.size() != knn)
    indices.set_size(knn);
  if (dists.size() != knn)
    dists.set_size(knn);

  flann::Matrix<float> query_mat ((float*)data_ptr, rows, cols); // Wrap query data
  flann::Matrix<int> index_mat(&indices[0], 1, knn);          // Wrap index vector

  if (m_dist_type == FLANN_DistType_L2) {

    // Allocate a temporary float distance matrix
    flann::Matrix<float> dists_mat (new float[dists.size()], 1, dists.size());
    int num_found = cast_index_ptr_L2_f(this->m_index_ptr)->
      knnSearch(query_mat, index_mat, dists_mat, knn, flann::SearchParams(128));

    // Copy the distances to the output matrix, then delete the temporary matrix
    for (size_t i=0; i<dists.size(); ++i)
      dists[i] = static_cast<double>(dists_mat[0][i]);
    delete[] dists_mat.ptr();
    return num_found;
  }

  vw_throw(IOErr() << "FLANNTree: Illegal distance type passed in.");
  return 0;
}

// Logic to be used for float and double types
template <class FLOAT_T>
void construct_index_aux(void* data_ptr, size_t num_features, size_t cols,
                         std::string const& flann_method, 
                         FLANN_DistType const& dist_type,
                         void* & index_ptr) {
  
  // When having a lot of features, and the user did not set a method, use the
  // kdtree method, which is 6x faster than kmeans but not deterministic. 
  int thresh = 25000; 
  std::string local_flann_method = flann_method; 
  if (local_flann_method == "auto" || local_flann_method == "") {
    if (num_features > thresh) {
      local_flann_method = "kdtree";
      {
        // Put a lock
        Mutex::WriteLock write_lock(flann_mutex);
        if (!flann_message_printed) {
          vw_out(vw::WarningMessage)            
             << "FLANNTree: Using the kdtree method for matching features with FLANN "
             << "as their number is more than " << thresh 
             << ". See the option --flann-method in the documentation.\n";
          flann_message_printed = true;
        }
      }
    } else {
     local_flann_method = "kmeans";
    }
  }
  
  if (index_ptr != NULL)
    vw_throw(IOErr() << "FLANNTree: Void ptr is not null, this is unexpected.");
  // Flann manual:
  // http://www.cs.ubc.ca/research/flann/uploads/FLANN/flann_manual-1.8.4.pdf.
  // Note: with Flann 1.9.2, the kdtree algorithm is not deterministic. It is
  // however about 6 times faster than the kmeans algorithm, which is
  // deterministic.  
  // Parameters for kdtree  
  const int NUM_TREES = 4;
  // Parameters for kmeans
  int branching = 32;
  int iterations = 11;
  flann::flann_centers_init_t centers_init = flann::FLANN_CENTERS_GONZALES;
  FLOAT_T cb_index = 0.2;
  switch (dist_type) {
  case FLANN_DistType_L2:
    if (local_flann_method == "kmeans") {
      // Slow, but deterministic
      index_ptr = new flann::Index<flann::L2<FLOAT_T>>(
          flann::Matrix<FLOAT_T>((FLOAT_T*)data_ptr, num_features, cols),
          flann::KMeansIndexParams(branching, iterations, centers_init, cb_index),
          flann::L2<FLOAT_T>());
    } else if (local_flann_method == "kdtree") {
      // Fast, but not deterministic
      index_ptr = new flann::Index<flann::L2<FLOAT_T>>(
          flann::Matrix<FLOAT_T>((FLOAT_T*)data_ptr, num_features, cols),
          flann::KDTreeIndexParams(NUM_TREES), // faster, but not deterministic
          flann::L2<FLOAT_T>());
    } else {
      vw_throw(IOErr() << "FLANNTree: Illegal FLANN method passed in.");
    }
    cast_index_ptr_L2_f(index_ptr)->buildIndex();
    return;
    
  case FLANN_DistType_Hamming:
    // Note that we use neither a kdtree or kmeans method
    index_ptr = new flann::Index<flann::Hamming<unsigned char>>
      (flann::Matrix<unsigned char>((unsigned char*)data_ptr, num_features, cols),
        //flann::LshIndexParams(), // Bad performance on small IP data sets
        flann::HierarchicalClusteringIndexParams(),
        flann::Hamming<unsigned char>());
    cast_index_ptr_HAMM_u(index_ptr)->buildIndex();
    return;
    
  default:
    vw_throw(IOErr() << "FLANNTree: Illegal distance type passed in.");
  }; // end switch
}

template <>
void FLANNTree<float>::construct_index(void* data_ptr, size_t num_features, size_t cols) {
  construct_index_aux<float>(data_ptr, num_features, cols, m_flann_method, m_dist_type, m_index_ptr);
}

template <>
FLANNTree<float>::~FLANNTree() {
  switch (m_dist_type) {
    case FLANN_DistType_L2:
      delete cast_index_ptr_L2_f(this->m_index_ptr); 
      return;
    default: 
      return;
  }; // end switch
}

template <>
size_t FLANNTree<float>::size1() const { 
  switch (m_dist_type) {
    case FLANN_DistType_L2:
      return cast_index_ptr_L2_f  (this->m_index_ptr)->size();
    default: 
      return 0;
  }; // end switch
}

template <>
size_t FLANNTree<float>::size2() const { 
  switch(m_dist_type) {
    case FLANN_DistType_L2:
      return cast_index_ptr_L2_f(this->m_index_ptr)->veclen();
    default: 
      return 0;
  }; // end switch
}

// All the same code duplicated but with doubles instead of floats
// TODO(oalexan1): Can this duplication be avoided?
template <>
size_t FLANNTree<double>::knn_search_help(void* data_ptr, size_t rows, size_t cols,
                                      Vector<int>& indices,
                                      Vector<double>& dists,
                                      size_t knn) {
  // Constrain the number of results that we can return to the number of loaded objects
  size_t maxNumReturns = m_features_cast.rows();
  if (knn > maxNumReturns)
    knn = maxNumReturns;

  if (indices.size() != knn)
    indices.set_size(knn);
  if (dists.size() != knn)
    dists.set_size(knn);

  flann::Matrix<double> query_mat ((double*)data_ptr, rows, cols); // Wrap query data
  flann::Matrix<int> index_mat(&indices[0], 1, knn);          // Wrap index vector
  if (m_dist_type == FLANN_DistType_L2) {

    flann::Matrix<double  > dists_mat (&dists[0], 1, dists.size()); // Wrap distance vector
    int num_found = cast_index_ptr_L2_d(this->m_index_ptr)->knnSearch
      (query_mat, index_mat, dists_mat, knn, flann::SearchParams(128));
    return num_found;
  }
  vw_throw(IOErr() << "FLANNTree: Illegal distance type passed in.");
  return 0;
}

// This code is duplicated for floats and doubles. 
template <>
void FLANNTree<double>::construct_index(void* data_ptr, size_t rows, size_t cols) {
  construct_index_aux<double>(data_ptr, rows, cols, m_flann_method, m_dist_type, m_index_ptr);
}

template <>
FLANNTree<double>::~FLANNTree() {
  switch(m_dist_type) {
    case FLANN_DistType_L2:
      delete cast_index_ptr_L2_d(this->m_index_ptr); return;
    default: 
      return;
  }; // end switch
}

template <>
size_t FLANNTree<double>::size1() const { 
  switch(m_dist_type) {
    case FLANN_DistType_L2:
      return cast_index_ptr_L2_d(this->m_index_ptr)->size();
    default: 
      return 0;
  }; // end switch
}

template <>
size_t FLANNTree<double>::size2() const { 
  switch(m_dist_type) {
    case FLANN_DistType_L2:
      return cast_index_ptr_L2_d(this->m_index_ptr)->veclen();
    default: 
      return 0;
  }; // end switch
}

//=============================================================================
// This is mostly the same code but uses Hamming distance instead of L2 distance

template <>
size_t FLANNTree<unsigned char>::knn_search_help(void* data_ptr, 
                                                size_t rows, size_t cols,
                                                Vector<int>& indices, 
                                                Vector<double>& dists,
                                                size_t knn) {
  // Constrain the number of results that we can return to the number of loaded objects
  size_t maxNumReturns = m_features_cast.rows();
  if (knn > maxNumReturns)
    knn = maxNumReturns;

  if (indices.size() != knn)
    indices.set_size(knn);
  if (dists.size() != knn)
    dists.set_size(knn);

  flann::Matrix<unsigned char> query_mat((unsigned char*)data_ptr, rows, cols); 
  flann::Matrix<int> index_mat(&indices[0], 1, knn); // Wrap index vector
  if (m_dist_type == FLANN_DistType_Hamming) {
    // Allocate a temporary uint32 distance matrix
    flann::Matrix<unsigned int> dists_mat (new unsigned int[knn], 1, knn);
    flann::SearchParams params;
    params.checks = 256; // Search more leaves
    params.cores  =   1; // Use only one core
    int num_found = cast_index_ptr_HAMM_u(this->m_index_ptr)->
      knnSearch(query_mat, index_mat, dists_mat, knn, params);
    // Copy the distances to the output matrix, then delete the temporary matrix
    for (size_t i=0; i<dists.size(); ++i)
      dists[i] = static_cast<double>(dists_mat[0][i]);
    delete[] dists_mat.ptr();
    return num_found;
  }

  vw_throw(IOErr() << "FLANNTree: Illegal distance type passed in.");
  return 0;
}

// TODO(oalexan1): Integrate this into construct_index_aux

template <>
void FLANNTree<unsigned char>::construct_index(void* data_ptr, size_t rows, size_t cols) {
  if (m_index_ptr != NULL)
    vw_throw(IOErr() << "FLANNTree: Void ptr is not null, this is unexpected.");
  switch(m_dist_type) {
    case FLANN_DistType_Hamming:
      m_index_ptr = new flann::Index<flann::Hamming<unsigned char>>
        (flann::Matrix<unsigned char>((unsigned char*)data_ptr, rows, cols),
          //flann::LshIndexParams(), // Bad performance on small IP data sets
          flann::HierarchicalClusteringIndexParams(),
          flann::Hamming<unsigned char>());
      cast_index_ptr_HAMM_u(this->m_index_ptr)->buildIndex();
      return;
  default:
    vw_throw(IOErr() << "FLANNTree: Illegal distance type passed in.");
  }; // end switch
}

template <>
FLANNTree<unsigned char>::~FLANNTree() {
  switch(m_dist_type) {
    case FLANN_DistType_Hamming: 
      delete cast_index_ptr_HAMM_u(this->m_index_ptr); 
      return;
    default: 
      return;
  }; // end switch
}

template <>
size_t FLANNTree<unsigned char>::size1() const { 
  switch(m_dist_type) {
    case FLANN_DistType_Hamming:
      return cast_index_ptr_HAMM_u(this->m_index_ptr)->size();
    default: 
      return 0;
  }; // end switch
}

template <>
size_t FLANNTree<unsigned char>::size2() const { 
  switch(m_dist_type) {
    case FLANN_DistType_Hamming:
      return cast_index_ptr_HAMM_u(this->m_index_ptr)->veclen();
    default: 
      return 0;
  }; // end switch
}

}}
