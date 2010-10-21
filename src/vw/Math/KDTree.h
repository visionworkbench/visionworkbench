// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


/// KD Tree
///
/// An implementation of a kd-tree, after "Multidimensional
/// Binary Search Trees Used for Associative Searching," Bentley 1975.
/// Implementation of m nearest neighbors follows "An Algorithm for finding best
/// Matches in Logarithmic Time," Friedman, Bently & Finkel
///
/// Tree only supports searches on records of numeric data
///
/// Methods
/// -- Building a KD tree (optimized, balanced, or randomly built,
///    depending on the the choice of partitioner and discriminator selector).
///   -- KD m nearest neighbors
///   -- Insertion of a single record into the tree
///   -- size of tree
/// -- Max and Min depth of tree (requires testing to insure that
///    depths are consistent regardless of whether build_tree or insert
///    was used.)
///   -- Region/Constrained Search (return all records within a region)
///      Constraint functors are likely to place constraints on records' keys,
///      but since the constraint functor has access to the record object, the
///      functor can be designed to access members other than a record's
///      iterable keys.
///
///
///  (TODO: list methods as they are implemented)
///
/// Desired Functionality to Implement:
///   -Deletion of the root
///   -Deletion of a random node
///   -Rebalance tree
///
///
/// KD-Tree Implementation Details:
///
/// Some terminology:
/// A file is a container of records. A record is a container of k keys.
///
/// Each vertex in the tree has associated properties, which can be accessed via
/// BGL property maps. The properties are:
///    -a record of k keys,
///    -a discriminator index,
///    -identifiers for its LO and HI subtrees, and
///    -hi and low ranges.
///
///
/// The record is the k-dimensional object. It must be iterable.
///
/// The discriminator identifies which of the vertex's k keys is used
/// to divide the subtrees. For some vertex P, let j =
/// discriminator(P). Then, any vertex Q in the LO subtree satisfies
/// Key_j(Q) < Key_j(P), and any vertex R in the HI subtree satisfies
/// Key_j(r) >= Key_j(Q). Note "=" for the HI subtree; this allows for
/// duplicate keys (multisets), but sacrifices any guarantees that the
/// tree be balanced.
///
/// Lo and Hi identifiers - bgl vertex descriptors 'loson' and 'hison'
/// which are the roots of the low and high subtrees below a vertex.
///
/// Lo and Hi ranges - Each vertex represents a partitioning of the
/// k-dimensional space along the dimension specified by the
/// discriminator. For any vertex P in the tree, all values in the
/// subtree rooted at P fall within the hyper-rectangle bounded by
/// lorange(P) and hirange(p)
///
///
///
///   The root's range is -infinity to infinity on each dimension, represented
///   as LO_Range = record of -inf, HI_Range = record of +inf.
///
///     L: [-inf, -inf, ..., -inf]
///     H: [+inf, +inf, ..., +inf]      root, disc = 1
///     Record:[0, 1, ..., -3]          | |
///                                    |   |
///                                   |     |
///       L: [-inf, -inf,...,-inf]   A        B  L:   [-inf, 1, ..., -inf]
///       H: (+inf, 1, ... ,+inf)                     H: [+inf, +inf, .., +inf]
///
///     By convention, the left has "<" and the rigt has ">="
///
///     Range variables should use the typedef range_t.
///
#ifndef __VW_MATH_KDTREE_H__
#define __VW_MATH_KDTREE_H__

#include <vw/Core/FundamentalTypes.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

#if (BOOST_VERSION >= 104000)
#include <boost/property_map/property_map.hpp>
#else
#include <boost/property_map.hpp>
#endif

// Std C++
#include <iostream>
#include <vector>
#include <iterator>
#include <queue>

// for INT_MAX
#include <climits>

//for printing graph
//#include <boost/graph/breadth_first_search.hpp>

//Create custom property tags. See section 3.6 in the BGL user guide
namespace boost{
  enum vertex_record_t{vertex_record = 511}; //A unique number
  BOOST_INSTALL_PROPERTY(vertex, record);

  enum vertex_discriminator_t{vertex_discriminator = 512};
  BOOST_INSTALL_PROPERTY(vertex, discriminator);

  enum vertex_LOSON_t{vertex_LOSON = 513};
  BOOST_INSTALL_PROPERTY(vertex, LOSON);

  enum vertex_HISON_t{vertex_HISON = 514};
  BOOST_INSTALL_PROPERTY(vertex, HISON);

  enum vertex_LORANGE_t{vertex_LORANGE = 515};
  BOOST_INSTALL_PROPERTY(vertex, LORANGE);

  enum vertex_HIRANGE_t{vertex_HIRANGE = 516};
  BOOST_INSTALL_PROPERTY(vertex, HIRANGE);
}

namespace vw {
namespace math {
  template<typename ForwardIterator>
  void print_file(ForwardIterator, ForwardIterator);
  template<typename ForwardIterator>
  void print_record(ForwardIterator, ForwardIterator);

  // Function object for comparing records by their discriminator keys
  // ContainerT must be have a const_iterator
  template<typename ContainerT>
  class DiscriminatorCompare {

    ContainerT m_pivot;
    unsigned m_discriminator;
    typename ContainerT::const_iterator pivot_iter;

  public:
    DiscriminatorCompare(ContainerT pivot, unsigned discriminator)
      : m_discriminator(discriminator){
      set_pivot(pivot);
    }

    DiscriminatorCompare(unsigned discriminator)
      : m_discriminator(discriminator){
      set_invalid_pivot();
    }

    // Unary operator. x[discriminator] < pivot[discriminator]
    bool operator()(const ContainerT& x) const {
      typename ContainerT::const_iterator input_iter = x.begin();
      std::advance(input_iter, m_discriminator);
      return *input_iter < *pivot_iter;
    }

    // Binary operator. x[discriminator] < y[discriminator]
    bool operator()(const ContainerT& x, const ContainerT& y) const {
      typename ContainerT::const_iterator x_iter = x.begin();
      typename ContainerT::const_iterator y_iter = y.begin();
      std::advance(x_iter, m_discriminator);
      std::advance(y_iter, m_discriminator);
      return *x_iter < *y_iter;
    }

    void set_pivot(ContainerT new_pivot){
      m_pivot = new_pivot;
      pivot_iter = m_pivot.begin();
      std::advance(pivot_iter, m_discriminator);
    }

    void set_invalid_pivot() {
      m_pivot = ContainerT();
      pivot_iter = m_pivot.end();
    }

    void set_discriminator(unsigned new_disc){m_discriminator = new_disc;}
  };


  /////////////////// Partitioners ////////////////////////////////////////

  // Partitions a file into two roughly equal parts.
  // Sorts the records in a file according to the key selected by
  // discriminator. The 'median' (not necessarily a true median if there are
  // duplicate values) is used to partition the file.
  //
  // Returns location of the partition, such that all elements less than
  // the partition are between beg and partition, and all elements
  // greater than or equal are in partition to end
  //
  // TODO: sort and partition may be overkill. Try to implement this with
  // partial_sort.
  class MedianPartitioner
  {
  public:
    template<typename RandomAccessIterT>
    RandomAccessIterT operator() (RandomAccessIterT beg, RandomAccessIterT end, unsigned discriminator) const
    {
      typedef typename std::iterator_traits<RandomAccessIterT>::value_type record_type;

      DiscriminatorCompare<record_type> comparator(discriminator);
      std::sort(beg, end, comparator); //Could use partial_sort() to save time?

      //Select the middle point as a pivot
      unsigned num_records = distance(beg, end);
      unsigned offset = (num_records / 2);
      RandomAccessIterT pivot_position = RandomAccessIterT(beg);
      std::advance(pivot_position, offset);

      //Partition the file around pivot
      comparator.set_pivot(*pivot_position);
      return std::partition(beg, end, comparator);
    }
  };

  // A thread-safer random number generator
  class KDRandom{
  public:
    ssize_t operator()(ssize_t max) const {
      double tmp;
      tmp = static_cast<double>(rand())
        / static_cast<double>(RAND_MAX);
      return static_cast<ssize_t>(tmp*max);
    }
  };


  //Input is partitioned around a randomly selected pivot record
  class RandPartitioner {
  public:
    template <typename RandomAccessIterT>
    RandomAccessIterT operator() (RandomAccessIterT beg, RandomAccessIterT end, unsigned disc) const
    {
      typedef typename std::iterator_traits<RandomAccessIterT>::value_type record_type;

      unsigned num_records = std::distance(beg, end);
      unsigned offset = rand() % num_records;
      RandomAccessIterT patition_position = beg;
      std::advance(patition_position, offset);
      DiscriminatorCompare<record_type> comparator( *patition_position, disc);
      std::sort(beg, end, comparator);
      return std::partition(beg, end, comparator);
    }
  };

  //////////////// Discriminator Selectors //////////////////////

  //Selects the discriminator as the dimension of maximum variance
  class VarianceDiscSelector {
  public:
    template<typename RandomAccessIterT>
    unsigned operator() (RandomAccessIterT file_beg, RandomAccessIterT file_end, unsigned /*disc*/) const {
      typedef typename std::iterator_traits<RandomAccessIterT>::value_type record_t;
      typedef typename record_t::iterator record_iter_t;

      //the number of keys in a record
      unsigned k = std::distance( (*file_beg).begin(), (*file_beg).end() );

      std::vector<double> mean;
      std::vector<double> variance;
      variance.assign(k,0);
      mean.assign(k,0);

      //Calculate Mean.
      int num_points = 0;
      for(RandomAccessIterT temp = file_beg; temp != file_end; ++temp){
        ++num_points;
        //mean +=record
        std::transform( (*temp).begin(), (*temp).end(), mean.begin(), mean.begin(), std::plus<double>() );
      }
      if(num_points != 0){
        for(unsigned j = 0; j<k; ++j){
          mean[j] /= num_points;
        }
      }

      //Calculate the variance
      std::vector<double>::iterator v_iter, m_iter;

      for( ; file_beg != file_end; ++file_beg){
        v_iter = variance.begin();
        m_iter = mean.begin();
        for(record_iter_t r_iter = (*file_beg).begin(); r_iter != (*file_beg).end(); ++r_iter){
          *v_iter += std::pow( (*r_iter - *m_iter), 2);
          ++v_iter;
          ++m_iter;
        }
      }

      //identify dimension of max variance
      double max_variance = -ScalarTypeLimits<double>::highest();
      unsigned index_of_max_variance = 0;
      for(unsigned j = 0; j<k; ++j){
        if(variance[j] > max_variance){
          max_variance = variance[j];
          index_of_max_variance = j;
        }
      }
      return index_of_max_variance;
    }
  };


  //Selects the discriminator as the dimension with largest (max-min) value
  class MaxDiffDiscSelector {
  public:
    template<typename RandomAccessIterT>
    unsigned operator() (RandomAccessIterT file_beg, RandomAccessIterT file_end, unsigned /*unused_argument*/) const {
      typedef typename std::iterator_traits<RandomAccessIterT>::value_type record_t;
      typedef typename record_t::iterator record_iter_t;

      unsigned k = std::distance( (*file_beg).begin(), (*file_beg).end() );

      //TODO: mins, maxes could be templated on the key type, and use
      //vw::ScalarTypeLimits<T>::highest()
      std::vector<double> mins;
      std::vector<double> maxes;
      mins.assign(k,ScalarTypeLimits<double>::highest());
      maxes.assign(k,-ScalarTypeLimits<double>::highest());

      std::vector<double> diffs;
      diffs.resize(k);

      unsigned i;
      for( ; file_beg != file_end; ++file_beg){
        i = 0;
        for(record_iter_t r_iter = (*file_beg).begin(); r_iter != (*file_beg).end(); ++r_iter){
          if(*r_iter < mins[i]) mins[i] = *r_iter;
          if(*r_iter > maxes[i]) maxes[i] = *r_iter;
          ++i;
        }
      }

      for(i=0; i<k; ++i){
        //TODO: make sure no wrap-around errors are possible
        diffs[i] = (maxes[i] - mins[i]);
      }

      //Identify index of max element
      double max_diff = -ScalarTypeLimits<double>::highest();
      unsigned dimension = 0;
      for(i=0; i<k; ++i){
        if(diffs[i] > max_diff){
          max_diff = diffs[i];
          dimension = i;
        }
      }
      return dimension;
    }
  };


  //NEXTDISC(d) = d+1 mod k (Bently's Discriminator selector)
  class ModuloDiscSelector {
    unsigned m_k;
  public:
    //Constructor
    ModuloDiscSelector(unsigned k)
      : m_k(k) {}
    template <typename ForwardIterT>
    unsigned operator() (ForwardIterT /*beg*/, ForwardIterT /*end*/, int disc) const {
      return (disc + 1) % m_k;
    }
  };


  ////////////////////////// Record Constraint Functors ////////////////////////
  // These must implement the unary operator(), and the binary
  // function domains overlap both functions must be passed objects with
  // iterators
  struct NullRecordConstraintKD {
    template<typename T>
    bool operator()(T /*record*/) const {
      return true;
    }

    template<typename T>
    bool domains_overlap(const T& /*lowRange*/, const T& /*highRange*/ ) const {
      return true;
    }
  };

  // When applied to a nearest neighbors search, restricts the result to
  // the m nearest records within the specified region
  template<typename RangeT>
  class RegionRecordConstraintKD {

    //Bounds of the constraint region:
    RangeT lowRange_;
    RangeT highRange_;

  public:
    //Constructor
    RegionRecordConstraintKD(RangeT lowRange, RangeT highRange)
      : lowRange_(lowRange), highRange_(highRange){}

    bool operator()(RangeT const& record) const {
      //check that all keys fall between lowRange_ and hiRange_
      typename RangeT::const_iterator low = lowRange_.begin();
      typename RangeT::const_iterator high = highRange_.begin();
      for(typename RangeT::const_iterator record_beg = record.begin(); record_beg != record.end(); ++record_beg)
        {
          if( *record_beg < *low )
            return false;
          if( *record_beg > *high )
            return false;
          ++low;
          ++high;
        }
      return true;
    }
    //Test if any part of a vertex's domain falls withiin the constrained region:
    bool domains_overlap(const RangeT &lowRange, const RangeT &highRange ) const {
      //check if either lowRange or highRange are within the constraint region
      //...odd syntax for calling this object's own operator()...
      bool lowWithin = operator()(lowRange);
      bool highWithin = operator()(highRange);
      if(lowWithin || highWithin)
        return true;
      //check if the constraint region is entirely within the input region
      RegionRecordConstraintKD inputRegionConstraint(lowRange, highRange);
      bool constraintLowWithin = inputRegionConstraint(lowRange_);
      bool constraintHighWithin = inputRegionConstraint(highRange_);
      if(constraintLowWithin && constraintHighWithin)
        return true;
      //else
      return false;
    }

  };

  //////////////////// Distance Metrics ////////////////////////////////////
  class SafeEuclideanDistanceMetric{
  public:
    //TODO: should a and b be required to be the same type of iterator?
    // If they were required to be of the same type, this could return
    //std::iterator_traits<IterT>::value_type
    //WARNING: may return max value for doubles
    /*
      template<typename AtypeT, typename BtypeT>
      double operator()(AtypeT a, BtypeT b){
      return operator()(a.begin(), a.end(), b.begin());
      }
    */
    template <typename IterA, typename IterB>
    double operator()(IterA a_first, IterA a_last, IterB b_first) const {
      double distance_sq = 0;

      for( ; a_first != a_last; ++a_first, ++b_first){
        distance_sq += pow( (*a_first - *b_first), 2);
      }
      return sqrt(distance_sq);
    }
  };


  /////////////////////////////////////////////////////////////////////////
  /// KD Tree
  //
  //  template <typename RandomAccessIterT>
  template <class FileT>
  class KDTree{

    typedef typename FileT::value_type record_t;
    typedef typename record_t::const_iterator record_iter_t;
    typedef typename std::iterator_traits<record_iter_t>::value_type key_t;

    typedef typename std::vector<key_t> range_t; //range_t should be selected to provide the operator[]

    typedef boost::graph_traits< boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS > >::vertex_descriptor Vertex;
    //      typedef boost::graph_traits< adjacency_list<vecS, vecS, undirectedS > >::edge_descriptor Edge;

    //The underlying graph structure:
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                  boost::property<boost::vertex_discriminator_t, unsigned,
                                    boost::property<boost::vertex_record_t, record_t,
                                             boost::property<boost::vertex_LOSON_t, Vertex,
                                                      boost::property<boost::vertex_HISON_t, Vertex,
                                                               boost::property<boost::vertex_LORANGE_t, range_t,
                                                                        boost::property<boost::vertex_HIRANGE_t, range_t > > > > > > > BGL_KDTree;

    //////////////////  --- KD Members --- /////////////////////////

    key_t m_POSITIVE_INFINITY;
    key_t m_NEGATIVE_INFINITY;
    BGL_KDTree m_kdTree;
    Vertex m_NIL;
    Vertex m_root;
    unsigned m_k; //dimensionality of kd tree
    unsigned m_max_depth;
    unsigned m_min_depth;

    class VertexDistanceComparator;

    //The m closest records encountered so far during a nearest neighbor search
    //The double is the distance of the Vertex from the query
    std::priority_queue< std::pair<Vertex, key_t>, std::vector<std::pair<Vertex, key_t> >, VertexDistanceComparator > m_priority_queue;

    // Property maps:
    // given a vertex descriptor V, its properties can be accessed like:
    // unsigned disc = m_discriminator_map[V];
    typename boost::property_map<BGL_KDTree, boost::vertex_record_t>::type m_record_map;
    typename boost::property_map<BGL_KDTree, boost::vertex_discriminator_t>::type m_discriminator_map;
    typename boost::property_map<BGL_KDTree, boost::vertex_LOSON_t>::type m_LOSON_map;
    typename boost::property_map<BGL_KDTree, boost::vertex_HISON_t>::type m_HISON_map;
    typename boost::property_map<BGL_KDTree, boost::vertex_LORANGE_t>::type m_LORANGE_map;
    typename boost::property_map<BGL_KDTree, boost::vertex_HIRANGE_t>::type m_HIRANGE_map;

  public:

    ///////////////  --- KD Constructors ---  //////////////////////

    KDTree(unsigned k) : m_k(k), m_max_depth(0), m_min_depth(INT_MAX){
      //initialize
      initialize_property_maps();
      initialize_infinity();

      //build tree:
      m_NIL = add_vertex(m_kdTree);
      m_root = m_NIL;
    }

    KDTree(unsigned k, FileT const& file) : m_k(k), m_max_depth(0), m_min_depth(INT_MAX)
    {
      //initialize
      initialize_property_maps();
      initialize_infinity();

      //build tree
      initialize_tree(file, ModuloDiscSelector(m_k), MedianPartitioner());
    }

    //Specify functors to select discriminator and partition
    template<typename DiscSelector, typename Partitioner>
    KDTree(unsigned k, FileT const& file,
           DiscSelector discSelector, Partitioner partitioner) : m_k(k), m_max_depth(0), m_min_depth(INT_MAX)
    {
      //initialize
      initialize_property_maps();
      initialize_infinity();

      //build tree
      initialize_tree(file, discSelector, partitioner);
    }

  private:

    //Capture the redundant initialization in the different constructors:
    void initialize_property_maps(){
      m_record_map = get(boost::vertex_record, m_kdTree);
      m_discriminator_map = get(boost::vertex_discriminator, m_kdTree);
      m_LOSON_map = get(boost::vertex_LOSON, m_kdTree);
      m_HISON_map = get(boost::vertex_HISON, m_kdTree);
      m_LORANGE_map = get(boost::vertex_LORANGE, m_kdTree);
      m_HIRANGE_map = get(boost::vertex_HIRANGE, m_kdTree);
    }
    void initialize_infinity(){
      m_POSITIVE_INFINITY = vw::ScalarTypeLimits<key_t>::highest();
      m_NEGATIVE_INFINITY = vw::ScalarTypeLimits<key_t>::lowest();
    }
    template<typename DiscSelector, typename Partitioner>
    void initialize_tree(FileT const& file,
                         DiscSelector /*discselector*/, Partitioner /*partitioner*/)
    {
      range_t lo_range(m_k, m_NEGATIVE_INFINITY);
      range_t hi_range(m_k, m_POSITIVE_INFINITY);

      m_NIL = add_vertex(m_kdTree);

      // TODO: This copy is required as the paritioners will run an
      // in-place sort on sections of the input.
      std::vector<record_t> temp_file(std::distance(file.begin(), file.end()));
      std::copy(file.begin(), file.end(), temp_file.begin());

      m_root = build_tree(temp_file.begin(), temp_file.end(), lo_range, hi_range,
                          ModuloDiscSelector(m_k), MedianPartitioner(), -1, 0);
    }


    //A private class provides the needed binary predicate for the queue:
    //With this, priority_queue.top() always returns the pair of greatest distance.
    class VertexDistanceComparator{
    public:
      template<typename T>
      bool operator()(const std::pair<Vertex,T>& pairA,
                      const std::pair<Vertex,T>& pairB) const{
        return (pairA.second < pairB.second);
      }
    };

  public:

    //////////////////   --- KD Public Methods ---  /////////////////////////

    //  Insert one record into an existing k-d tree. Does not guarantee
    //  a balanced tree. Allowing duplicate records produces a multiset, where all records
    //  in the left subtree are < root, and all records in the right subtree are >= the root.
    void insert(record_t r){
      kd_insert(r, m_root, 0);
    }

    /*
    //TODO:
    void clear(){
    //remove all vertices
    //clear all state
    }
    */


    /// M_NEAREST_NEIGHBORS
    //
    // Fills nearest_records with (up to) m records which are nearest to query.
    // the return int is the number of records returned, which will be m unless
    // the tree contains fewer than m records, in which case the number of records in
    // the tree is returned.
    //
    // Query must be a container with k keys.
    //TODO: the return type of nearest_records should not depend on ContainerT
    template <typename ContainerT>
    unsigned m_nearest_neighbors(ContainerT const& query,
                                 std::vector<record_t>& nearest_records,
                                 unsigned m = 1)
    {
      return m_nearest_neighbors(query, nearest_records, m,
                                 NullRecordConstraintKD(),
                                 SafeEuclideanDistanceMetric());
    }

    template <typename ContainerT, typename RecordConstraintT, typename DistanceMetricT>
    unsigned m_nearest_neighbors(ContainerT const& query,
                                 std::vector<record_t>& nearest_records,
                                 unsigned m = 1,
                                 RecordConstraintT recordConstraint = NullRecordConstraintKD(),
                                 DistanceMetricT distanceMetric = SafeEuclideanDistanceMetric())
    {
      assert( m_k == (unsigned) std::distance(query.begin(), query.end()) );
      nearest_records.clear();

      //Initialize the priority queue with m sentinel NIL pairs
      std::pair<Vertex, key_t> nil_pair(m_NIL, m_POSITIVE_INFINITY);
      while (!m_priority_queue.empty()){
        m_priority_queue.pop();
      }
      for(unsigned i = 0; i < m; ++i){
        m_priority_queue.push(nil_pair);
      }

      //convert input into a range_t object, which is more convenient for
      //internal operations (provides the operator[])
      range_t _query(m_k);
      std::copy(query.begin(), query.end(), _query.begin());

      nearest_neighbors(m_root, _query, m, recordConstraint, distanceMetric);

      //Now parse m_priority_queue into nearest_records
      unsigned num_records_found = 0;
      std::pair<Vertex, key_t> temp_pair;

      while (!m_priority_queue.empty()){

        temp_pair = m_priority_queue.top();
        if (temp_pair.first == m_NIL){
          m_priority_queue.pop();
        }else{
          //add record of neighbor node
          ++num_records_found;
          nearest_records.push_back(m_record_map[temp_pair.first]);
          m_priority_queue.pop();
        }
      }
      // The priority queue is ordered in decreasing distance from the query,
      // so nearest_records will need to be reversed after it is filled
      reverse (nearest_records.begin(), nearest_records.end());

      return num_records_found;
    }

    //The number of vertices in the graph, less one for the NIL vertex
    unsigned size(){return num_vertices(m_kdTree)-1;}
    //The root vertex is at depth 0.
    //TODO: ...but an empty tree also has depth 0... hmmm...
    unsigned max_depth(){return m_max_depth;}
    unsigned min_depth(){return m_min_depth;}

  private:

    ///////////////// --- KD Private Methods --- ////////////////////

    /// BUILD_TREE
    //
    // Given all the records, building the k-d tree is a matter of
    // selecting the proper discriminator and a pivot which evenly partitions
    // the input file. Returns the root of the k-d tree.
    //
    // DiscSelector:a functor which must support
    //  unsigned operator() (IterT beg, IterT end, unsigned disc)
    //  and return a discriminator between 0 and k-1, inclusive.
    //
    // Partitioner: a functor which must support
    //  IterT operator() (IterT beg, IterT end, unsigned disc)
    //  and return the position of an element to be used for partitioning
    //
    template <typename IterT, typename DiscSelector, typename Partitioner>
    Vertex build_tree(IterT file_beg, IterT file_end,
                      range_t lo_range, range_t hi_range,
                      DiscSelector discSelector, Partitioner partitioner,
                      int previous_disc, unsigned depth)
    {
      // End recursion?
      if (file_beg == file_end){
        if(depth < m_min_depth){
          m_min_depth = depth;
        }
        return m_NIL;
      }

      if (depth > m_max_depth){
        m_max_depth = depth;
      }

      //Choose discriminator and partition
      unsigned disc = discSelector(file_beg, file_end, previous_disc);
      IterT partition = partitioner(file_beg, file_end, disc);

      Vertex P = add_vertex(m_kdTree);
      m_record_map[P] = *partition;
      m_discriminator_map[P] = disc;

      m_LORANGE_map[P] = lo_range;
      m_HIRANGE_map[P] = hi_range;

      //Create high range for LOSON:
      range_t LOSON_hi_range = hi_range;
      set_range_value(LOSON_hi_range.begin(), (*partition).begin(), disc);
      //Create low range for HISON
      range_t HISON_lo_range = lo_range;
      set_range_value(HISON_lo_range.begin(), (*partition).begin(), disc);

      //The subtrees below P should not include P
      IterT one_past_partition = IterT(partition);
      std::advance(one_past_partition,1);

      Vertex lo = build_tree(file_beg, partition, lo_range, LOSON_hi_range, discSelector, partitioner, disc, ++depth);
      Vertex hi = build_tree(one_past_partition, file_end, HISON_lo_range, hi_range, discSelector, partitioner, disc, ++depth);
      m_LOSON_map[P] = lo;
      m_HISON_map[P] = hi;

      //TODO:
      //Also insert bgl edges. This increases memory usage, so perhaps
      //should be commented out after debugging?
      /*
        Edge ed;
        bool inserted;
        boost::tie(ed, inserted) = add_edge(P,lo,m_kdTree);
        boost::tie(ed, inserted) = add_edge(P,hi,m_kdTree);
      */

      //Print out P, LOSON, HISON, discriminator, ranges
      //print__vertex(P, lo, hi);
      return P;
    }


    // Insert attempts to find a vertex with a nil child. The edge
    // (both the explicit BGL edge and the implicit edge defined by LOSON[P]
    // or HISON[P] ) from parent to the nil child are removed. A new vertex
    // is created for the new record, and it is inserted as a child to repace
    // P's previous nil child.
    // The new vertex's discriminator is determined using via modulo
    //
    // TODO: test updates for m_min_depth and m_max_depth
    void kd_insert(record_t r, Vertex P, unsigned depth){
      if(P == m_NIL){
        //this must be an empty tree
        Vertex new_vertex = add_vertex(m_kdTree);
        m_min_depth = 0;
        m_max_depth = 0;

        //initialize properties for new vertex
        //range_t lo_range(m_k, -ScalarTypeLimits<double>::highest());
        //range_t hi_range(m_k, ScalarTypeLimits<double>::highest());
        range_t lo_range(m_k, m_NEGATIVE_INFINITY);
        range_t hi_range(m_k, m_POSITIVE_INFINITY);

        m_record_map[new_vertex] = r;
        m_discriminator_map[new_vertex] = 0;
        m_LORANGE_map[new_vertex] = lo_range;
        m_HIRANGE_map[new_vertex] = hi_range;
        m_LOSON_map[new_vertex] = m_NIL;
        m_HISON_map[new_vertex] = m_NIL;

        //BGL Edges
        /*
          Edge ed;
          bool inserted;
          boost::tie(ed, inserted) = add_edge(new_vertex, m_NIL,m_kdTree);
          boost::tie(ed, inserted) = add_edge(new_vertex, m_NIL,m_kdTree);
        */

        m_root = new_vertex;
        return;
      }

      unsigned disc_P = m_discriminator_map[P];

      if( on_low_side(r, P) ){
        if( m_LOSON_map[P] == m_NIL ){
          //Insert r as P's loson:
          //remove_edge(P, m_NIL, m_kdTree); //explicit BGL edge
          Vertex new_vertex = add_vertex(m_kdTree);
          m_LOSON_map[P] = new_vertex;  //implicit edge

          //initialize properties for new vertex
          record_t record_P = m_record_map[P];
          range_t new_hirange = m_HIRANGE_map[P];
          set_range_value(new_hirange.begin(), record_P.begin(), disc_P);

          m_record_map[new_vertex] = r;
          m_discriminator_map[new_vertex] = ( (disc_P + 1) % m_k);
          m_HIRANGE_map[new_vertex] = new_hirange;
          m_LORANGE_map[new_vertex] = m_LORANGE_map[P];
          m_LOSON_map[new_vertex] = m_NIL;
          m_HISON_map[new_vertex] = m_NIL;

          //BGL Edges
          /*
            Edge ed;
            bool inserted;
            boost::tie(ed, inserted) = add_edge(new_vertex, m_NIL, m_kdTree);
            boost::tie(ed, inserted) = add_edge(new_vertex, m_NIL, m_kdTree);
          */

          //P is at depth 'depth', so P's son is at depth 'depth+1'
          if(depth+1 < m_min_depth)
            m_min_depth = depth+1;

          if (depth+1 > m_max_depth)
            m_max_depth = depth+1;

          return;
        }
        //otherwise continue the recursion
        kd_insert(r, m_LOSON_map[P], ++depth);
      }
      else{ //on high side of partition
        if(m_HISON_map[P] == m_NIL){
          // Insert r as P's hison
          //remove_edge(P, m_NIL, m_kdTree);
          Vertex new_vertex = add_vertex(m_kdTree);
          m_HISON_map[P] = new_vertex;

          // Initialize properties for new vertex
          record_t record_P = m_record_map[P];
          range_t new_lorange = m_LORANGE_map[P];
          set_range_value(new_lorange.begin(), record_P.begin(), disc_P);

          m_record_map[new_vertex] = r;
          m_discriminator_map[new_vertex] = ( (disc_P + 1) % m_k);
          m_LORANGE_map[new_vertex] = new_lorange;
          m_HIRANGE_map[new_vertex] = m_HIRANGE_map[P];
          m_LOSON_map[new_vertex] = m_NIL;
          m_HISON_map[new_vertex] = m_NIL;

          //BGL Edges
          /*
            Edge ed;
            bool inserted;
            boost::tie(ed, inserted) = add_edge(new_vertex, m_NIL, m_kdTree);
            boost::tie(ed, inserted) = add_edge(new_vertex, m_NIL, m_kdTree);
          */

          //P is at depth 'depth', so P's son is at depth 'depth+1'
          if(depth+1 < m_min_depth)
            m_min_depth = depth+1;

          if (depth+1 > m_max_depth)
            m_max_depth = depth+1;

          return;
        }
        //otherwise continue the recursion
        kd_insert(r, m_HISON_map[P], ++depth);
      }
    }


    /// Nearest Neighbors
    //
    // Finds the m nearest records to the query, and stores them in
    // m_priority_queue.
    //
    // The return value 0 indicates 'return' from recursion, whereas a
    // return value of 1 indicates the search is complete
    // TODO: Consider using std::make_pair in place of boost's tie function
    //
    // Constrained searches would be faster if children's ranges were checked
    // for overlap with the constraint region: only children whose ranges
    // might contain records satisfying the constraint would be explored.
    // This could be done by only searching a child if the child's lorange or
    // the child's hirange satisfies the constraint.
    template<typename RecordConstraint, typename DistanceMetric>
    int nearest_neighbors(Vertex const& N, range_t const& query, unsigned m,
                          RecordConstraint const& recordConstraint,
                          DistanceMetric const& distanceMetric ) {
      if (N == m_NIL)
        return 0;

      range_t& N_lorange = get(m_LORANGE_map, N);
      range_t& N_hirange = get(m_HIRANGE_map, N);

      //If no part of N's domain satisfies the constraint, don't search
      if(!recordConstraint.domains_overlap(N_lorange,  N_hirange)){
        return 0;
      }
      record_t& N_record = get(m_record_map, N);
      unsigned disc = get(m_discriminator_map, N);
      Vertex& loson = get(m_LOSON_map, N);
      Vertex& hison = get(m_HISON_map, N);

      //Calculate distance from query to N's record.
      //double distance_to_N = kd_euclidean_distance(query.begin(), query.end(), N_record.begin());
      double distance_to_N = distanceMetric(query.begin(), query.end(), N_record.begin());

      //Update the list of m best.
      Vertex dummy = m_NIL;
      double distance_to_mth_best = m_POSITIVE_INFINITY;
      if (!m_priority_queue.empty()){
        boost::tie(dummy, distance_to_mth_best) = m_priority_queue.top();
      } else {
        std::cout<<"Priority Queue is empty!?\n";
      }

      if ((distance_to_N < distance_to_mth_best) && recordConstraint(N_record) ){
        //if (distance_to_N < distance_to_mth_best){
        std::pair<Vertex, double> temp_pair(N, distance_to_N);
        m_priority_queue.push(temp_pair);

        while(m_priority_queue.size() > m){
          m_priority_queue.pop();
        }
        boost::tie(dummy, distance_to_mth_best) = m_priority_queue.top();
      }

      // Now it is necesary to search the subtree on the same side of
      // the partition as the query. When recursion returns, if it is
      // possible that the other subtree could contain nearer
      // neighbors, it is necessary to search that subtree too.
      record_iter_t partition_key_position = N_record.begin();
      std::advance(partition_key_position, disc);
      double partition_key = *partition_key_position;
      double query_key = query[disc];

      int done = 0;
      if (query_key < partition_key)
        {
          //query is on low side of partition at vertex N
          done = nearest_neighbors(loson, query, m, recordConstraint, distanceMetric);
          if (done == 1)
            return 1;
          //m nearest may have changed from searching loson
          if (!m_priority_queue.empty())
            boost::tie(dummy, distance_to_mth_best) = m_priority_queue.top();
          N_lorange[disc] = partition_key;
          if (bounds_overlap_ball(query, N_lorange, N_hirange, distance_to_mth_best))
            {
              //std::cout<<"Hison's Bounds Overlap Ball, check hi side too.\n";
              done = nearest_neighbors(hison, query, m, recordConstraint, distanceMetric);
              if(done == 1)
                return 1;
            }
        }else{
        //query is on the high side of the partition
        done = nearest_neighbors(hison, query, m, recordConstraint, distanceMetric);
        if (done == 1)
          return 1;
        //m nearest may have changed from searching hison
        if (!m_priority_queue.empty())
          boost::tie(dummy, distance_to_mth_best) = m_priority_queue.top();
        N_hirange[disc] = partition_key;
        if (bounds_overlap_ball(query, N_lorange, N_hirange, distance_to_mth_best))
          {
            //std::cout<<"LOson's's Bounds Overlap Ball, check lo side too.\n";
            done = nearest_neighbors(loson, query, m, recordConstraint, distanceMetric);
            if (done == 1)
              return 1;
          }
      }
      //After checking children, the ball containing the m best may be entirely in N's range:
      boost::tie(dummy, distance_to_mth_best) = m_priority_queue.top();
      if(ball_within_bounds(query, N_lorange, N_hirange, distance_to_mth_best) ){
        return 1;
      }
      return 0;
    }



    // The following two functions refer to a "ball" centered at the query.
    // The radius of the ball is the distance from the query to the mth
    // nearest record examined so far in the search. The ball, then,
    // represents a region where possibly closer records could be found. If a
    // region overlaps this ball, that means part of the ball is within that
    // region's bounds, and the region must be
    // examined. (bounds_overlap_ball) If a ball fits entirely within a
    // region that has been examined, the search is
    // finished. (ball_within_bounds)


    /// Bounds Overlap Ball
    //
    // Checks if ANY part of the ball lies within the given bounds.
    // If the distance from the query to the bounded region is greater
    // than the radius then there is no overlap. This is implemented by
    // checking if d^2 = d_1^2 + d_2^2 + ... + d_k^2 > r^2
    //
    // If the radius is infinite, than it is considered to overlap any
    // bounds, even if they are also infinite.
    inline bool bounds_overlap_ball(range_t const& query, range_t const& lorange, range_t const& hirange, double radius) const {

      //Now we can assume that the query has no infinities in it.
      double radius_squared = pow(radius, 2);
      double d_squared = 0;
      typename range_t::const_iterator iter = query.begin();
      typename range_t::const_iterator lorange_iter = lorange.begin();
      typename range_t::const_iterator hirange_iter = hirange.begin();
      for( ; iter != query.end(); ++iter, ++lorange_iter, ++hirange_iter) {

        if(*iter < *lorange_iter) {
          d_squared += pow( *iter - *lorange_iter, 2);
          if (d_squared > radius_squared)
            return false;

        } else if(*iter > *hirange_iter) {
          d_squared += pow( *iter - *hirange_iter, 2);
          if (d_squared > radius_squared)
            return false;
        }
        //else if *lorange_iter < *iter < *hirange_iter then within
        //bounds, and the ith coordiante distance is 0
      }
      //d_squared < radius_squared, so no overlap
      return true;
    }

    /// Ball Within Bounds
    //
    // Checks if a sphere of given radius fits entirely within the bounds
    // by checking if the 'coordinate distance' from the query to the
    // bounded region is less than the radius of the ball.
    //
    // TODO: all uses of ScalarTypeLimits<double>::highest() should be replaced with a C implementation of infinity!
    bool ball_within_bounds(range_t const& query, range_t const& lorange, range_t const& hirange, double radius){

      //If any coordinate distance is less less than radius, ball
      //escapes bounds, so return false if a bound is +/- infinity,
      //the coordinate distance is infinity
      typename range_t::const_iterator iter = query.begin();
      typename range_t::const_iterator lorange_iter = lorange.begin();
      typename range_t::const_iterator hirange_iter = hirange.begin();
      for( ; iter != query.end(); ++iter, ++lorange_iter, ++hirange_iter) {
        // DELETE?
        //         if ((fabs(*lorange_iter) == m_POSITIVE_INFINITY) || (fabs(*hirange_iter) == m_POSITIVE_INFINITY))
        //           return false;
        if (fabs( *iter - *lorange_iter ) <= radius)
          return false;
        if (fabs( *iter - *hirange_iter ) <= radius)
          return false;
      }
      return true;
    }

    //checks if a record is on the low side of the vertex's partition.
    //That is, if d is V's discriminator, record[d] < V[d]
    //V probably shouldn't be NIL...
    template<typename QueryT>
    bool on_low_side(QueryT r, Vertex V)
    {
      unsigned disc = m_discriminator_map[V];
      record_t partition = m_record_map[V];
      record_iter_t partition_iter = partition.begin();
      typename QueryT::iterator r_iter = r.begin();
      std::advance(partition_iter, disc);
      std::advance(r_iter, disc);
      return *r_iter < *partition_iter;
    }

    //Performs A[k] = B[k]
    template<typename IterA, typename IterB>
    void set_range_value(IterA iterA, IterB iterB, unsigned k){
      std::advance(iterA, k);
      std::advance(iterB, k);
      *iterA = *iterB;
    }

    //Given vertex descriptors for a vertex Q, its loson P and hison R,
    //output the vertices' records, and the discriminator used
    void print_vertex(Vertex P, Vertex Q, Vertex R){
      record_t P_record = get(m_record_map, P);
      record_t Q_record = get(m_record_map, Q);
      record_t R_record = get(m_record_map, R);
      range_t P_lorange = get(m_LORANGE_map, P);
      range_t P_hirange = get(m_HIRANGE_map, P);

      std::cout<< "\tRecord[P] : ";
      print_record(P_record.begin(), P_record.end());

      std::cout<< "\tRecord[loson] : ";
      print_record(Q_record.begin(), Q_record.end());

      std::cout<< "\tRecord[hison] : ";
      print_record(R_record.begin(), R_record.end());

      std::cout<<"\tDiscrimnator[P] : " <<m_discriminator_map[P]<<"\n";
      std::cout<<"Lo-Range of Subtree rooted at P: ";
      print_record(P_lorange.begin(), P_lorange.end());
      std::cout<<"Hi-Range of Subtree rooted at P: ";
      print_record(P_hirange.begin(), P_hirange.end());
      std::cout<<"\n";
    }
  }; //end KDTree class


  //Display all of the records in a file.
  template<typename ForwardIterator>
  void print_file (ForwardIterator begin_records, ForwardIterator end_records)
  {
    for( ; begin_records != end_records; ++begin_records){
      print_record((*begin_records).begin(), (*begin_records).end());
    }
  }

  //Display the contents of a container (e.g. a record)
  template <typename ForwardIterator>
  void print_record(ForwardIterator first, ForwardIterator last){
    //Identify "infinity"
    typedef typename std::iterator_traits<ForwardIterator>::value_type key_t;
    key_t positive_infinity = vw::ScalarTypeLimits<key_t>::highest();
    key_t negative_infinity = vw::ScalarTypeLimits<key_t>::lowest();
    std::cout<<"[ ";
    for( ; first !=last; ++first){
      if(*first == positive_infinity){
        std::cout<<"+inf ";
      }else if(*first == negative_infinity){
        std::cout<<"-inf ";
      }else{
        std::cout<<*first<< " ";
      }
    }
    std::cout<<" ]\n";
  }

}} // namespace vw::math
#endif // __VW_MATH_KDTREE_H__


