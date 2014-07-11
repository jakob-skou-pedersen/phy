/******************************************************************** 
 * Copyright (C) 2008-2014 Jakob Skou Pedersen - All Rights Reserved.
 *
 * See README_license.txt for license agreement.
 *******************************************************************/
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE phyTestRunner
#include <boost/test/unit_test.hpp>

#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

using boost::unit_test::test_suite;

// BOOST_AUTO_TEST_CASE(my_test_function_1) 
// {
//   BOOST_CHECK( 2 == 1 );
// }
// 
// BOOST_AUTO_TEST_CASE(my_test_function_2) 
// {
//   BOOST_CHECK( 2 == 2 );
// }
// 
// typedef boost::mpl::list<int,long,unsigned char> test_types;
// BOOST_AUTO_TEST_CASE_TEMPLATE( my_test, T, test_types )
// {
//   BOOST_CHECK_EQUAL( sizeof(T), (unsigned)4 );
// }
// 
// 
// // test suites: http://www.boost.org/doc/libs/1_36_0/libs/test/doc/html/utf/user-guide/test-organization/auto-test-suite.html
// BOOST_AUTO_TEST_SUITE( test_suite1 )
// 
//   BOOST_AUTO_TEST_CASE( test_case1 )
// {
//   BOOST_WARN( sizeof(int) < 4 );
// }
// 
// BOOST_AUTO_TEST_CASE( test_case2 )
// {
//   BOOST_REQUIRE_EQUAL( 1, 2 );
//   BOOST_FAIL( "Should never reach this line" );
// }
// 
// BOOST_AUTO_TEST_SUITE_END();
// 
// BOOST_AUTO_TEST_SUITE( test_suite2 )
// 
// BOOST_AUTO_TEST_CASE( test_case3 )
// {
//   BOOST_CHECK( true );
// }
// 
// BOOST_AUTO_TEST_CASE( test_case4 )
// {
//   BOOST_CHECK( false );
// }
// 
// BOOST_AUTO_TEST_SUITE_END()
