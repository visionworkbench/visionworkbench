/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ErrorPrinter.h>

int main() {
 return CxxTest::ErrorPrinter().run();
}
#include "TestGPU_ImageStatistics.h"

static TestImageStatistics suite_TestImageStatistics;

static CxxTest::List Tests_TestImageStatistics = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestImageStatistics( "TestGPU_ImageStatistics.h", 35, "TestImageStatistics", suite_TestImageStatistics, Tests_TestImageStatistics );

static class TestDescription_TestImageStatistics_test_min_channel_value : public CxxTest::RealTestDescription {
public:
 TestDescription_TestImageStatistics_test_min_channel_value() : CxxTest::RealTestDescription( Tests_TestImageStatistics, suiteDescription_TestImageStatistics, 39, "test_min_channel_value" ) {}
 void runTest() { suite_TestImageStatistics.test_min_channel_value(); }
} testDescription_TestImageStatistics_test_min_channel_value;

static class TestDescription_TestImageStatistics_test_max_channel_value : public CxxTest::RealTestDescription {
public:
 TestDescription_TestImageStatistics_test_max_channel_value() : CxxTest::RealTestDescription( Tests_TestImageStatistics, suiteDescription_TestImageStatistics, 50, "test_max_channel_value" ) {}
 void runTest() { suite_TestImageStatistics.test_max_channel_value(); }
} testDescription_TestImageStatistics_test_max_channel_value;

static class TestDescription_TestImageStatistics_test_min_max_channel_values : public CxxTest::RealTestDescription {
public:
 TestDescription_TestImageStatistics_test_min_max_channel_values() : CxxTest::RealTestDescription( Tests_TestImageStatistics, suiteDescription_TestImageStatistics, 61, "test_min_max_channel_values" ) {}
 void runTest() { suite_TestImageStatistics.test_min_max_channel_values(); }
} testDescription_TestImageStatistics_test_min_max_channel_values;

static class TestDescription_TestImageStatistics_test_mean_channel_value : public CxxTest::RealTestDescription {
public:
 TestDescription_TestImageStatistics_test_mean_channel_value() : CxxTest::RealTestDescription( Tests_TestImageStatistics, suiteDescription_TestImageStatistics, 74, "test_mean_channel_value" ) {}
 void runTest() { suite_TestImageStatistics.test_mean_channel_value(); }
} testDescription_TestImageStatistics_test_mean_channel_value;

#include <cxxtest/Root.cpp>
