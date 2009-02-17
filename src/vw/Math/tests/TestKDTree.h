// __BEGIN_LICENSE__
// Copyright (C) 2006, 2007 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <cxxtest/TestSuite.h>
#include <vw/Math/KDTree.h>

using namespace vw;
using namespace vw::math;

class TestKDTree : public CxxTest::TestSuite
{
  public:

  void test_discriminator_compare()
    {
      std::vector<int> recordA, recordB, recordC;

      recordA.push_back(1);
      recordA.push_back(2);
      recordA.push_back(3);

      recordB.push_back(-2);
      recordB.push_back(2);
      recordB.push_back(2);

      recordC.push_back(10);
      recordC.push_back(12);
      recordC.push_back(5);


      DiscriminatorCompare<std::vector<int> > comparator(recordA, 0);
      //test unary operator() input[disc] < pivot[disc]
      TS_ASSERT(comparator(recordB));
      TS_ASSERT(!comparator(recordC));
      comparator.set_discriminator(1);
      TS_ASSERT(!comparator(recordB));

      //test binary operator()
      TS_ASSERT(comparator(recordA, recordC));
      TS_ASSERT(!comparator(recordA, recordA));

      comparator.set_pivot(recordC);
      assert(comparator(recordA));
    }

  void test_median_partitioner(){
    MedianPartitioner mp;

    std::vector< std::vector<int> > file(6);
    typedef std::vector< std::vector<int> >::iterator file_iter_t;

    file[0].push_back(1);
    file[0].push_back(2);
    file[0].push_back(3);

    file[1].push_back(7);
    file[1].push_back(10);
    file[1].push_back(-1);

    file[2].push_back(8);
    file[2].push_back(18);
    file[2].push_back(0);

    file[3].push_back(-22);
    file[3].push_back(18);
    file[3].push_back(9);

    file[4].push_back(0);
    file[4].push_back(0);
    file[4].push_back(0);

    file[5].push_back(-2);
    file[5].push_back(-2);
    file[5].push_back(-2);

    //partition with discriminator = 0
    unsigned disc = 0;
    file_iter_t partition_position = mp(file.begin(), file.end(), disc);
    std::vector<int> partition = *partition_position;
    TS_ASSERT_EQUALS(partition[0], 1);
    TS_ASSERT_EQUALS(partition[1], 2);
    TS_ASSERT_EQUALS(partition[2], 3);

    for(file_iter_t iter = file.begin(); iter != partition_position; ++iter){
      TS_ASSERT( (*iter)[disc] < (*partition_position)[disc] );
    }

    for(file_iter_t iter = partition_position; iter !=file.end(); ++iter){
      TS_ASSERT( (*iter)[disc] >= (*partition_position)[disc] );
    }


    //partition with discriminator = 1
    disc = 1;
    partition_position = mp(file.begin(), file.end(), disc);
    partition = *partition_position;
    TS_ASSERT_EQUALS(partition[0], 7);
    TS_ASSERT_EQUALS(partition[1], 10);
    TS_ASSERT_EQUALS(partition[2], -1);

    for(file_iter_t iter = file.begin(); iter != partition_position; ++iter){
      TS_ASSERT( (*iter)[disc] < (*partition_position)[disc] );
    }

    for(file_iter_t iter = partition_position; iter !=file.end(); ++iter){
      TS_ASSERT( (*iter)[disc] >= (*partition_position)[disc] );
    }


    //partition with discriminator = 2
    disc = 2;
    partition_position = mp(file.begin(), file.end(), disc);
    partition = *partition_position;
    TS_ASSERT_EQUALS(partition[disc], 0);

    for(file_iter_t iter = file.begin(); iter != partition_position; ++iter){
      TS_ASSERT( (*iter)[disc] < (*partition_position)[disc] );
    }

    for(file_iter_t iter = partition_position; iter !=file.end(); ++iter){
      TS_ASSERT( (*iter)[disc] >= (*partition_position)[disc] );
    }


  }

  void test_rand_partitioner(){
    //The rand partitioner could be checked by verifying that
    //all records on the low side of the partition satisfy
    //low_record[disc] < partition[disc] and all records on hi
    //side satisfy hi_record[disc] >= partition[disc]
    std::vector< std::vector<double> > file(8);
    typedef std::vector<std::vector<double> >::iterator file_iter_t;

    file[0].push_back(0.1);
    file[0].push_back(2.2);
    file[0].push_back(10);

    file[1].push_back(0.2);
    file[1].push_back(-2.8);
    file[1].push_back(-10);

    file[2].push_back(-0.1);
    file[2].push_back(1.4);
    file[2].push_back(9);

    file[3].push_back(0.23);
    file[3].push_back(3.1);
    file[3].push_back(13);

    file[4].push_back(0.2);
    file[4].push_back(2);
    file[4].push_back(11.5);

    file[5].push_back(0);
    file[5].push_back(-2.9);
    file[5].push_back(-17);

    file[6].push_back(0.1);
    file[6].push_back(1.4);
    file[6].push_back(9);

    file[7].push_back(-0.23);
    file[7].push_back(2.1);
    file[7].push_back(-15);

    RandPartitioner rp;
    unsigned disc = 0;
    file_iter_t partition = rp(file.begin(), file.end(), disc);
    for(file_iter_t iter = file.begin(); iter != partition; ++iter){
      TS_ASSERT( (*iter)[disc] < (*partition)[disc] );
    }

    for(file_iter_t iter = partition; iter !=file.end(); ++iter){
      TS_ASSERT( (*iter)[disc] >= (*partition)[disc] );
    }


  }

  void test_variance_disc_selector(){
    std::vector< std::vector<double> > fileA(8);

    fileA[0].push_back(0.1);
    fileA[0].push_back(2.2);
    fileA[0].push_back(10);

    fileA[1].push_back(0.2);
    fileA[1].push_back(-2.8);
    fileA[1].push_back(-10);

    fileA[2].push_back(-0.1);
    fileA[2].push_back(1.4);
    fileA[2].push_back(9);

    fileA[3].push_back(0.23);
    fileA[3].push_back(3.1);
    fileA[3].push_back(13);

    fileA[4].push_back(0.2);
    fileA[4].push_back(2);
    fileA[4].push_back(11.5);

    fileA[5].push_back(0);
    fileA[5].push_back(-2.9);
    fileA[5].push_back(-17);

    fileA[6].push_back(0.1);
    fileA[6].push_back(1.4);
    fileA[6].push_back(9);

    fileA[7].push_back(-0.23);
    fileA[7].push_back(2.1);
    fileA[7].push_back(-15);

    VarianceDiscSelector vds;
    //The int argument does nothing
    unsigned dimension_of_max_variance = vds(fileA.begin(), fileA.end(), 0);
    TS_ASSERT_EQUALS(int(dimension_of_max_variance), 2);

    std::vector< std::vector<double> > fileB(3);

    fileB[0].push_back(0.1);
    fileB[0].push_back(222.2);
    fileB[0].push_back(10);
    fileB[0].push_back(10);

    fileB[1].push_back(0.2);
    fileB[1].push_back(-235.8);
    fileB[1].push_back(-10);
    fileB[1].push_back(10);


    fileB[2].push_back(-0.1);
    fileB[2].push_back(155.4);
    fileB[2].push_back(9);
    fileB[2].push_back(10);


    dimension_of_max_variance = vds(fileB.begin(), fileB.end(), 22);
    TS_ASSERT_EQUALS(int(dimension_of_max_variance), 1);
  }

  void test_modulo_disc_selector(){
    ModuloDiscSelector mds(4);

    std::vector< std::vector<double> > file(2);

    file[0].push_back(0.1);
    file[0].push_back(222.2);
    file[0].push_back(10);
    file[0].push_back(10);

    file[1].push_back(0.2);
    file[1].push_back(-235.8);
    file[1].push_back(-10);
    file[1].push_back(10);

    TS_ASSERT_EQUALS(int(mds(file.begin(), file.end(), -1)), 0);
    TS_ASSERT_EQUALS(int(mds(file.begin(), file.end(), 0)), 1);
    TS_ASSERT_EQUALS(int(mds(file.begin(), file.end(), 1)), 2);
    TS_ASSERT_EQUALS(int(mds(file.begin(), file.end(), 2)), 3);
    TS_ASSERT_EQUALS(int(mds(file.begin(), file.end(), 3)), 0);
  }

  void test_max_diff_disc_selector(){
    MaxDiffDiscSelector mdds;

    std::vector< std::vector<double> > fileB(3);

    fileB[0].push_back(0.1);
    fileB[0].push_back(222.2);
    fileB[0].push_back(10);
    fileB[0].push_back(10);

    fileB[1].push_back(0.2);
    fileB[1].push_back(-235.8);
    fileB[1].push_back(-10);
    fileB[1].push_back(10);


    fileB[2].push_back(-0.1);
    fileB[2].push_back(155.4);
    fileB[2].push_back(9);
    fileB[2].push_back(10);

    TS_ASSERT_EQUALS(int(mdds(fileB.begin(), fileB.end(), 7)), 1);
  }

  void test_basic_constructor(){
    typedef std::vector< std::vector<int> > file_t;
    KDTree<file_t> kd(5);
    TS_ASSERT_EQUALS(int(kd.size()), 0);

    //insert
    std::vector<int> new_record;

    new_record.push_back(4);
    new_record.push_back(5);
    new_record.push_back(6);
    new_record.push_back(5);
    new_record.push_back(5);

    kd.insert(new_record);
    TS_ASSERT_EQUALS(int(kd.size()), 1);
  }

  void test_constructor_default_functors(){

    std::vector< std::vector<int> > file(6);
    typedef std::vector< std::vector<int> > file_t;

    file[0].push_back(1);
    file[0].push_back(2);
    file[0].push_back(3);

    file[1].push_back(7);
    file[1].push_back(10);
    file[1].push_back(-1);

    file[2].push_back(8);
    file[2].push_back(18);
    file[2].push_back(0);

    file[3].push_back(-22);
    file[3].push_back(18);
    file[3].push_back(9);

    file[4].push_back(0);
    file[4].push_back(0);
    file[4].push_back(0);

    file[5].push_back(-2);
    file[5].push_back(-2);
    file[5].push_back(-2);

    KDTree<file_t> kd(3, file);
    TS_ASSERT_EQUALS(int(kd.size()), 6);

    //insert
    std::vector<int> new_record;
    new_record.push_back(4);
    new_record.push_back(5);
    new_record.push_back(6);

    kd.insert(new_record);
    TS_ASSERT_EQUALS(int(kd.size()), 7);
  }

  void test_search(){

    typedef std::vector< std::vector<double> > file_t;
    file_t file(6);

    file[0].push_back(1.0);
    file[0].push_back(2.2);
    file[0].push_back(3.3);

    file[1].push_back(7);
    file[1].push_back(10);
    file[1].push_back(-1);

    file[2].push_back(8);
    file[2].push_back(18);
    file[2].push_back(0);

    file[3].push_back(-22);
    file[3].push_back(18);
    file[3].push_back(9);

    file[4].push_back(0);
    file[4].push_back(0);
    file[4].push_back(0);

    file[5].push_back(-2);
    file[5].push_back(-2);
    file[5].push_back(-2);

    KDTree<file_t> kd_A(3, file);
    KDTree<file_t> kd_B(3, file, VarianceDiscSelector(), MedianPartitioner());
    KDTree<file_t> kd_C(3, file, VarianceDiscSelector(), RandPartitioner());
    KDTree<file_t> kd_D(3, file, ModuloDiscSelector(3), RandPartitioner());
    KDTree<file_t> kd_E(3, file, MaxDiffDiscSelector(), MedianPartitioner());


    std::vector<double> queryA, queryB, queryC, queryD;
    queryA.push_back(-0.5);
    queryA.push_back(-0.5);
    queryA.push_back(-0.5);

    queryB.push_back(-50);
    queryB.push_back(50);
    queryB.push_back(10);

    queryC.push_back(-3);
    queryC.push_back(-3);
    queryC.push_back(-3);

    queryD.push_back(20);
    queryD.push_back(20);
    queryD.push_back(20);

    std::vector<std::vector<double> > nearest_records;
    int num_records;

    ////////////////////////////////////////////////////////////////////////////////
    //    TS_TRACE("Search in tree A.");

    //    TS_TRACE("Find 4 nearest to query ");
    //    print_record(queryB.begin(), queryB.end());
    num_records = kd_A.m_nearest_neighbors(queryB, nearest_records, 4);
    TS_ASSERT_EQUALS(num_records, 4);

    //    print_file(nearest_records.begin(), nearest_records.end());

    TS_ASSERT_EQUALS(nearest_records[0][0], -22);
    TS_ASSERT_EQUALS(nearest_records[0][1], 18);
    TS_ASSERT_EQUALS(nearest_records[0][2], 9);

    //TS_TRACE(stringify("Size of kd_A: ") + stringify(kd_A.size()));
    TS_ASSERT_EQUALS(nearest_records[1][0], 8); //err
    TS_ASSERT_EQUALS(nearest_records[1][1], 18);//err
    TS_ASSERT_EQUALS(nearest_records[1][2], 0); //err

    //TS_TRACE("Find 2 nearest to query ")
    //print_record(queryA.begin(), queryA.end());
    num_records = kd_A.m_nearest_neighbors(queryA, nearest_records, 2);
    TS_ASSERT_EQUALS(num_records, 2);
    TS_ASSERT_EQUALS(nearest_records[0][0], 0);
    TS_ASSERT_EQUALS(nearest_records[0][1], 0);
    TS_ASSERT_EQUALS(nearest_records[0][2], 0);

    TS_ASSERT_EQUALS(nearest_records[1][0], -2);
    TS_ASSERT_EQUALS(nearest_records[1][1], -2);
    TS_ASSERT_EQUALS(nearest_records[1][2], -2);

    //searching for more records than are in the tree should return the
    //size of the tree
    //TS_TRACE("Find 14 nearest to query);
    //print_record(queryA.begin(), queryA.end());
    num_records = kd_A.m_nearest_neighbors(queryA, nearest_records, 14);
    TS_ASSERT_EQUALS(num_records, 6);
    TS_ASSERT_EQUALS(nearest_records[0][0], 0);
    TS_ASSERT_EQUALS(nearest_records[0][1], 0);
    TS_ASSERT_EQUALS(nearest_records[0][2], 0);

    TS_ASSERT_EQUALS(nearest_records[1][0], -2);
    TS_ASSERT_EQUALS(nearest_records[1][1], -2);
    TS_ASSERT_EQUALS(nearest_records[1][2], -2);

    //TS_TRACE("Find 1 nearest to query");
    //print_record(queryC.begin(), queryC.end());
    num_records = kd_A.m_nearest_neighbors(queryC, nearest_records);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], -2);
    TS_ASSERT_EQUALS(nearest_records[0][1], -2);
    TS_ASSERT_EQUALS(nearest_records[0][2], -2);

    /////////////////////////////////////////////////////////////////////////
    //TS_TRACE("Search in tree B.");

    //TS_TRACE("Find 2 nearest to query")
    //print_record(queryA.begin(), queryA.end());
    num_records = kd_B.m_nearest_neighbors(queryA, nearest_records, 2);
    TS_ASSERT_EQUALS(num_records, 2);
    TS_ASSERT_EQUALS(nearest_records[0][0], 0);
    TS_ASSERT_EQUALS(nearest_records[0][1], 0);
    TS_ASSERT_EQUALS(nearest_records[0][2], 0);

    TS_ASSERT_EQUALS(nearest_records[1][0], -2);
    TS_ASSERT_EQUALS(nearest_records[1][1], -2);
    TS_ASSERT_EQUALS(nearest_records[1][2], -2);

    //searching for more records than are in the tree should return the
    //size of the tree

    //TS_TRACE("Find 14 nearest to query");
    //print_record(queryA.begin(), queryA.end());
    num_records = kd_B.m_nearest_neighbors(queryA, nearest_records, 14);
    TS_ASSERT_EQUALS(num_records, 6);
    TS_ASSERT_EQUALS(nearest_records[0][0], 0);
    TS_ASSERT_EQUALS(nearest_records[0][1], 0);
    TS_ASSERT_EQUALS(nearest_records[0][2], 0);

    TS_ASSERT_EQUALS(nearest_records[1][0], -2);
    TS_ASSERT_EQUALS(nearest_records[1][1], -2);
    TS_ASSERT_EQUALS(nearest_records[1][2], -2);

    //TS_TRACE("Find 3 nearest to query");
    //print_record(queryB.begin(), queryB.end());
    num_records = kd_B.m_nearest_neighbors(queryB, nearest_records, 3);
    TS_ASSERT_EQUALS(num_records, 3);
    TS_ASSERT_EQUALS(nearest_records[0][0], -22);
    TS_ASSERT_EQUALS(nearest_records[0][1], 18);
    TS_ASSERT_EQUALS(nearest_records[0][2], 9);

    TS_ASSERT_EQUALS(nearest_records[1][0], 8); //err
    TS_ASSERT_EQUALS(nearest_records[1][1], 18); //err
    TS_ASSERT_EQUALS(nearest_records[1][2], 0); //err

    //TS_TRACE("Find 1 nearest to query");
    //print_record(queryC.begin(), queryC.end());
    num_records = kd_B.m_nearest_neighbors(queryC, nearest_records);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], -2);
    TS_ASSERT_EQUALS(nearest_records[0][1], -2);
    TS_ASSERT_EQUALS(nearest_records[0][2], -2);



    /////////////////////////////////////////////////////////////////////////
    //TS_TRACE("Search in tree C.");

    num_records = kd_C.m_nearest_neighbors(queryA, nearest_records, 2);
    TS_ASSERT_EQUALS(num_records, 2);
    TS_ASSERT_EQUALS(nearest_records[0][0], 0);
    TS_ASSERT_EQUALS(nearest_records[0][1], 0);
    TS_ASSERT_EQUALS(nearest_records[0][2], 0);

    TS_ASSERT_EQUALS(nearest_records[1][0], -2);
    TS_ASSERT_EQUALS(nearest_records[1][1], -2);
    TS_ASSERT_EQUALS(nearest_records[1][2], -2);

    //searching for more records than are in the tree should return the
    //size of the tree
    num_records = kd_C.m_nearest_neighbors(queryA, nearest_records, 14);
    TS_ASSERT_EQUALS(num_records, 6);
    TS_ASSERT_EQUALS(nearest_records[0][0], 0);
    TS_ASSERT_EQUALS(nearest_records[0][1], 0);
    TS_ASSERT_EQUALS(nearest_records[0][2], 0);

    TS_ASSERT_EQUALS(nearest_records[1][0], -2);
    TS_ASSERT_EQUALS(nearest_records[1][1], -2);
    TS_ASSERT_EQUALS(nearest_records[1][2], -2);

    num_records = kd_C.m_nearest_neighbors(queryB, nearest_records, 3);
    TS_ASSERT_EQUALS(num_records, 3);
    TS_ASSERT_EQUALS(nearest_records[0][0], -22);
    TS_ASSERT_EQUALS(nearest_records[0][1], 18);
    TS_ASSERT_EQUALS(nearest_records[0][2], 9);

    TS_ASSERT_EQUALS(nearest_records[1][0], 8);
    TS_ASSERT_EQUALS(nearest_records[1][1], 18);
    TS_ASSERT_EQUALS(nearest_records[1][2], 0);

    num_records = kd_C.m_nearest_neighbors(queryC, nearest_records);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], -2);
    TS_ASSERT_EQUALS(nearest_records[0][1], -2);
    TS_ASSERT_EQUALS(nearest_records[0][2], -2);


    /////////////////////////////////////////////////////////////////////
    //TS_TRACE("Search in tree D.");

    num_records = kd_D.m_nearest_neighbors(queryA, nearest_records, 2);
    TS_ASSERT_EQUALS(num_records, 2);
    TS_ASSERT_EQUALS(nearest_records[0][0], 0);
    TS_ASSERT_EQUALS(nearest_records[0][1], 0);
    TS_ASSERT_EQUALS(nearest_records[0][2], 0);

    TS_ASSERT_EQUALS(nearest_records[1][0], -2);
    TS_ASSERT_EQUALS(nearest_records[1][1], -2);
    TS_ASSERT_EQUALS(nearest_records[1][2], -2);

    //searching for more records than are in the tree should return the
    //size of the tree
    num_records = kd_D.m_nearest_neighbors(queryA, nearest_records, 14);
    TS_ASSERT_EQUALS(num_records, 6);
    TS_ASSERT_EQUALS(nearest_records[0][0], 0);
    TS_ASSERT_EQUALS(nearest_records[0][1], 0);
    TS_ASSERT_EQUALS(nearest_records[0][2], 0);

    TS_ASSERT_EQUALS(nearest_records[1][0], -2);
    TS_ASSERT_EQUALS(nearest_records[1][1], -2);
    TS_ASSERT_EQUALS(nearest_records[1][2], -2);

    num_records = kd_D.m_nearest_neighbors(queryB, nearest_records, 3);
    TS_ASSERT_EQUALS(num_records, 3);
    TS_ASSERT_EQUALS(nearest_records[0][0], -22);
    TS_ASSERT_EQUALS(nearest_records[0][1], 18);
    TS_ASSERT_EQUALS(nearest_records[0][2], 9);

    TS_ASSERT_EQUALS(nearest_records[1][0], 8);
    TS_ASSERT_EQUALS(nearest_records[1][1], 18);
    TS_ASSERT_EQUALS(nearest_records[1][2], 0);

    num_records = kd_D.m_nearest_neighbors(queryC, nearest_records);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], -2);
    TS_ASSERT_EQUALS(nearest_records[0][1], -2);
    TS_ASSERT_EQUALS(nearest_records[0][2], -2);


    ////////////////////////////////////////////////////////////////////
    //TS_TRACE("Search in tree E.");

    num_records = kd_E.m_nearest_neighbors(queryA, nearest_records, 2);
    TS_ASSERT_EQUALS(num_records, 2);
    TS_ASSERT_EQUALS(nearest_records[0][0], 0);
    TS_ASSERT_EQUALS(nearest_records[0][1], 0);
    TS_ASSERT_EQUALS(nearest_records[0][2], 0);

    TS_ASSERT_EQUALS(nearest_records[1][0], -2);
    TS_ASSERT_EQUALS(nearest_records[1][1], -2);
    TS_ASSERT_EQUALS(nearest_records[1][2], -2);

    //searching for more records than are in the tree should return the
    //size of the tree
    num_records = kd_E.m_nearest_neighbors(queryA, nearest_records, 14);
    TS_ASSERT_EQUALS(num_records, 6);
    TS_ASSERT_EQUALS(nearest_records[0][0], 0);
    TS_ASSERT_EQUALS(nearest_records[0][1], 0);
    TS_ASSERT_EQUALS(nearest_records[0][2], 0);

    TS_ASSERT_EQUALS(nearest_records[1][0], -2);
    TS_ASSERT_EQUALS(nearest_records[1][1], -2);
    TS_ASSERT_EQUALS(nearest_records[1][2], -2);

    num_records = kd_E.m_nearest_neighbors(queryB, nearest_records, 3);
    TS_ASSERT_EQUALS(num_records, 3);
    TS_ASSERT_EQUALS(nearest_records[0][0], -22);
    TS_ASSERT_EQUALS(nearest_records[0][1], 18);
    TS_ASSERT_EQUALS(nearest_records[0][2], 9);

    TS_ASSERT_EQUALS(nearest_records[1][0], 8); //err
    TS_ASSERT_EQUALS(nearest_records[1][1], 18); //err
    TS_ASSERT_EQUALS(nearest_records[1][2], 0); //err

    num_records = kd_E.m_nearest_neighbors(queryC, nearest_records);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], -2);
    TS_ASSERT_EQUALS(nearest_records[0][1], -2);
    TS_ASSERT_EQUALS(nearest_records[0][2], -2);
  }

  void test_searches_after_inserts(){

    typedef std::vector<double> record_t;
    std::vector<record_t> file(9);
    typedef std::vector< record_t > file_t;

    file[0].push_back(0);
    file[0].push_back(0);

    file[1].push_back(100);
    file[1].push_back(0);

    file[2].push_back(0);
    file[2].push_back(10);

    file[3].push_back(0);
    file[3].push_back(100);

    file[4].push_back(10);
    file[4].push_back(0);

    file[5].push_back(-10);
    file[5].push_back(-10);

    file[6].push_back(-10);
    file[6].push_back(-1000);

    file[7].push_back(10);
    file[7].push_back(0);

    file[8].push_back(10);
    file[8].push_back(10);

    //Some records to insert:
    record_t recordA, recordB, recordC;

    recordA.push_back(1000);
    recordA.push_back(1000);

    recordB.push_back(0);
    recordB.push_back(0);

    recordC.push_back(-10);
    recordC.push_back(-1001);

    //Some queries
    record_t queryA, queryB, queryC, queryD;

    queryA.push_back(0);
    queryA.push_back(0);

    queryB.push_back(3);
    queryB.push_back(3);

    queryC.push_back(-10);
    queryC.push_back(-1002);

    queryD.push_back(2000);
    queryD.push_back(2000);



    std::vector<record_t> nearest_records;
    int num_records;


    KDTree<file_t> kd_A(2);
    KDTree<file_t> kd_B(2, file);
    KDTree<file_t> kd_C(2, file, VarianceDiscSelector(), MedianPartitioner());
    KDTree<file_t> kd_D(2, file, VarianceDiscSelector(), RandPartitioner());
    KDTree<file_t> kd_E(2, file, ModuloDiscSelector(2), RandPartitioner());
    KDTree<file_t> kd_F(2, file, MaxDiffDiscSelector(), MedianPartitioner());

    //Attempt a search on an empty tree:
    num_records = kd_A.m_nearest_neighbors(queryB, nearest_records, 4);
    TS_ASSERT_EQUALS(num_records, 0);
    TS_ASSERT_EQUALS(kd_A.max_depth(), 0);

    //Insert vertices and the search
    kd_A.insert(recordA);
    num_records = kd_A.m_nearest_neighbors(queryB, nearest_records, 4);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(kd_A.max_depth(), 0);
    TS_ASSERT_EQUALS(nearest_records[0][0], recordA[0]);
    TS_ASSERT_EQUALS(nearest_records[0][1], recordA[1]);

    //insert another record.
    kd_A.insert(recordB);
    num_records = kd_A.m_nearest_neighbors(queryB, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(kd_A.max_depth(), 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], recordB[0]);
    TS_ASSERT_EQUALS(nearest_records[0][1], recordB[1]);

    //insert another record.
    kd_A.insert(recordC);
    num_records = kd_A.m_nearest_neighbors(queryB, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(kd_A.max_depth(), 2);
    TS_ASSERT_EQUALS(nearest_records[0][0], recordB[0]);
    TS_ASSERT_EQUALS(nearest_records[0][1], recordB[1]);

    //////////////////////   kd_B   /////////////////////////////////////

    num_records = kd_B.m_nearest_neighbors(queryC, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], -10);
    TS_ASSERT_EQUALS(nearest_records[0][1], -1000);

    kd_B.insert(recordC);
    num_records = kd_B.m_nearest_neighbors(queryC, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], recordC[0]);
    TS_ASSERT_EQUALS(nearest_records[0][1], recordC[1]);

    kd_B.insert(recordA);
    kd_B.insert(recordB);
    num_records = kd_B.m_nearest_neighbors(queryD, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], recordA[0]);
    TS_ASSERT_EQUALS(nearest_records[0][1], recordA[1]);


    //////////////////////   kd_C   /////////////////////////////////////

    num_records = kd_C.m_nearest_neighbors(queryC, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], -10);
    TS_ASSERT_EQUALS(nearest_records[0][1], -1000);

    kd_C.insert(recordC);
    num_records = kd_C.m_nearest_neighbors(queryC, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], recordC[0]);
    TS_ASSERT_EQUALS(nearest_records[0][1], recordC[1]);

    kd_C.insert(recordA);
    kd_C.insert(recordB);
    num_records = kd_C.m_nearest_neighbors(queryD, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], recordA[0]);
    TS_ASSERT_EQUALS(nearest_records[0][1], recordA[1]);


    //////////////////////   kd_D   /////////////////////////////////////

    num_records = kd_D.m_nearest_neighbors(queryC, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], -10);
    TS_ASSERT_EQUALS(nearest_records[0][1], -1000);

    kd_D.insert(recordC);
    num_records = kd_D.m_nearest_neighbors(queryC, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], recordC[0]);
    TS_ASSERT_EQUALS(nearest_records[0][1], recordC[1]);

    kd_D.insert(recordA);
    kd_D.insert(recordB);
    num_records = kd_D.m_nearest_neighbors(queryD, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], recordA[0]);
    TS_ASSERT_EQUALS(nearest_records[0][1], recordA[1]);


    //////////////////////   kd_E   /////////////////////////////////////

    num_records = kd_E.m_nearest_neighbors(queryC, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], -10);
    TS_ASSERT_EQUALS(nearest_records[0][1], -1000);

    kd_E.insert(recordC);
    num_records = kd_E.m_nearest_neighbors(queryC, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], recordC[0]);
    TS_ASSERT_EQUALS(nearest_records[0][1], recordC[1]);

    kd_E.insert(recordA);
    kd_E.insert(recordB);
    num_records = kd_E.m_nearest_neighbors(queryD, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], recordA[0]);
    TS_ASSERT_EQUALS(nearest_records[0][1], recordA[1]);


    //////////////////////   kd_F   /////////////////////////////////////

    num_records = kd_F.m_nearest_neighbors(queryC, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], -10);
    TS_ASSERT_EQUALS(nearest_records[0][1], -1000);

    kd_F.insert(recordC);
    num_records = kd_F.m_nearest_neighbors(queryC, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], recordC[0]);
    TS_ASSERT_EQUALS(nearest_records[0][1], recordC[1]);

    kd_F.insert(recordA);
    kd_F.insert(recordB);
    num_records = kd_F.m_nearest_neighbors(queryD, nearest_records, 1);
    TS_ASSERT_EQUALS(num_records, 1);
    TS_ASSERT_EQUALS(nearest_records[0][0], recordA[0]);
    TS_ASSERT_EQUALS(nearest_records[0][1], recordA[1]);

  }

  void test_region_record_constraint(){
    typedef std::vector<double> range_t;
    double infinity = ScalarTypeLimits<double>::highest();
    range_t lowRange, highRange;
    lowRange.push_back(0);
    lowRange.push_back(0);
    highRange.push_back(infinity);
    highRange.push_back(infinity);
    RegionRecordConstraintKD<range_t> regionRecordConstraintKD(lowRange, highRange);

    //some records
    range_t recordA, recordB, recordC;

    recordA.push_back(10);
    recordA.push_back(10);

    recordB.push_back(-10);
    recordB.push_back(10);

    recordC.push_back(0);
    recordC.push_back(0);


    //Record A is in the range specified
    TS_ASSERT(regionRecordConstraintKD(recordA));
    //Record B is not
    TS_ASSERT(!regionRecordConstraintKD(recordB));
    //Record C is on a boundary, which is accepted
    TS_ASSERT(regionRecordConstraintKD(recordC));

    //Test domains_overlap

    range_t recordD, recordE, recordF;

    recordD.push_back(-15);
    recordD.push_back(-5);

    recordE.push_back(-2);
    recordE.push_back(-2);

    recordF.push_back(-15);
    recordF.push_back(-5);


    //Part of the input region overlaps the constraint region:
    TS_ASSERT(regionRecordConstraintKD.domains_overlap(recordA, recordB));
    //No overlap between input region and constraint region
    TS_ASSERT(!regionRecordConstraintKD.domains_overlap(recordD, recordE));
    //Constraint region within input region, with equality or high bounds:
    TS_ASSERT(regionRecordConstraintKD.domains_overlap(recordD, highRange));
    //Input region entirely within constraint region
    TS_ASSERT(regionRecordConstraintKD.domains_overlap(recordA, recordA));

    //Constraint region entirely within input region:
    RegionRecordConstraintKD<range_t> anotherRegionRecordConstraintKD(recordC, recordA);
    TS_ASSERT(anotherRegionRecordConstraintKD.domains_overlap(recordD, highRange));

  }



  void test_safe_euclidean_distance_metric(){
    TS_WARN("Need to test safe euclidean distance metric");
  }

  void test_constrained_search(){
    typedef std::vector<double> range_t;
    double infinity = ScalarTypeLimits<double>::highest();
    range_t lowRange, highRange;
    lowRange.push_back(0);
    lowRange.push_back(0);
    highRange.push_back(infinity);
    highRange.push_back(infinity);
    //Constrain search to only find records in the first quadrant:
    RegionRecordConstraintKD<range_t> regionRecordConstraintKD(lowRange, highRange);


    //Build A KD tree
     std::vector< range_t > file(6);
    typedef std::vector< range_t > file_t;

    file[0].push_back(1.0);
    file[0].push_back(2.2);

    file[1].push_back(7);
    file[1].push_back(10);

    file[2].push_back(8);
    file[2].push_back(18);

    file[3].push_back(-22);
    file[3].push_back(18);

    file[4].push_back(0);
    file[4].push_back(0);

    file[5].push_back(-2);
    file[5].push_back(-2);

    KDTree<file_t> kd(2, file);

    //Some queries

    range_t queryA, queryB, queryC;

    queryA.push_back(-3);
    queryA.push_back(-3);

    queryB.push_back(-25);
    queryB.push_back(25);

    queryC.push_back(9);
    queryC.push_back(19);

    //Perform some constrained searches

    std::vector<range_t> nearest_records;
    int num_records;

    //Find nearest record, excluding those outside the search region:
    num_records = kd.m_nearest_neighbors(queryA, nearest_records, 1, regionRecordConstraintKD, SafeEuclideanDistanceMetric());
    TS_ASSERT_EQUALS( nearest_records[0][0], 0);
    TS_ASSERT_EQUALS( nearest_records[0][1], 0);

    //If any key fails the constraint, like the record (-22, 18), the record should be excluded:
    num_records = kd.m_nearest_neighbors(queryB, nearest_records, 1, regionRecordConstraintKD, SafeEuclideanDistanceMetric());
    TS_ASSERT_EQUALS( nearest_records[0][0], 8);
    TS_ASSERT_EQUALS( nearest_records[0][1], 18);

    //Searches within the constrained region should be unaffected:
    num_records = kd.m_nearest_neighbors(queryC, nearest_records, 1, regionRecordConstraintKD, SafeEuclideanDistanceMetric());
    TS_ASSERT_EQUALS( nearest_records[0][0], 8);
    TS_ASSERT_EQUALS( nearest_records[0][1], 18);

    //should find all records in the first quardant:
    num_records = kd.m_nearest_neighbors(file[5], nearest_records, 22, regionRecordConstraintKD, SafeEuclideanDistanceMetric());
    TS_ASSERT_EQUALS(num_records, 4);
    TS_ASSERT_EQUALS( nearest_records[0][0], 0);
    TS_ASSERT_EQUALS( nearest_records[0][1], 0);
    TS_ASSERT_EQUALS( nearest_records[1][0], 1);
    TS_ASSERT_EQUALS( nearest_records[1][1], 2.2);
    TS_ASSERT_EQUALS( nearest_records[2][0], 7);
    TS_ASSERT_EQUALS( nearest_records[2][1], 10);
    TS_ASSERT_EQUALS( nearest_records[3][0], 8);
    TS_ASSERT_EQUALS( nearest_records[3][1], 18);
  }

};
