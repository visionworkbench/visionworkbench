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
#include "TestGPU_Algorithms.h"

static TestAlgorithms suite_TestAlgorithms;

static CxxTest::List Tests_TestAlgorithms = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestAlgorithms( "TestGPU_Algorithms.h", 35, "TestAlgorithms", suite_TestAlgorithms, Tests_TestAlgorithms );

static class TestDescription_TestAlgorithms_test_image_algo_fill : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAlgorithms_test_image_algo_fill() : CxxTest::RealTestDescription( Tests_TestAlgorithms, suiteDescription_TestAlgorithms, 39, "test_image_algo_fill" ) {}
 void runTest() { suite_TestAlgorithms.test_image_algo_fill(); }
} testDescription_TestAlgorithms_test_image_algo_fill;

static class TestDescription_TestAlgorithms_test_image_algo_clamp : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAlgorithms_test_image_algo_clamp() : CxxTest::RealTestDescription( Tests_TestAlgorithms, suiteDescription_TestAlgorithms, 53, "test_image_algo_clamp" ) {}
 void runTest() { suite_TestAlgorithms.test_image_algo_clamp(); }
} testDescription_TestAlgorithms_test_image_algo_clamp;

static class TestDescription_TestAlgorithms_test_image_algo_normalize : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAlgorithms_test_image_algo_normalize() : CxxTest::RealTestDescription( Tests_TestAlgorithms, suiteDescription_TestAlgorithms, 85, "test_image_algo_normalize" ) {}
 void runTest() { suite_TestAlgorithms.test_image_algo_normalize(); }
} testDescription_TestAlgorithms_test_image_algo_normalize;

static class TestDescription_TestAlgorithms_test_image_algo_threshold : public CxxTest::RealTestDescription {
public:
 TestDescription_TestAlgorithms_test_image_algo_threshold() : CxxTest::RealTestDescription( Tests_TestAlgorithms, suiteDescription_TestAlgorithms, 117, "test_image_algo_threshold" ) {}
 void runTest() { suite_TestAlgorithms.test_image_algo_threshold(); }
} testDescription_TestAlgorithms_test_image_algo_threshold;

#include <cxxtest/Root.cpp>
