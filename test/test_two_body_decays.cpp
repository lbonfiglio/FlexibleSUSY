#include "two_body_decays.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_two_body_decays

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace flexiblesusy::two_body_decays;

BOOST_AUTO_TEST_CASE( test_kallen )
{
   BOOST_CHECK_EQUAL(kallen(0., 0., 0.), 0.);
   BOOST_CHECK_EQUAL(kallen(1.0, 1.0, 1.0), -3.);
   BOOST_CHECK_EQUAL(kallen(1.0, 2.0, 0.0), 1.);
   BOOST_CHECK_EQUAL(kallen(3.0, 0.0, 2.0), 1.);
}

BOOST_AUTO_TEST_CASE( test_scalar_to_fermion_fermion )
{
   BOOST_CHECK_EQUAL(scalar_to_fermion_fermion(100., 60., 50., 0.1, 0.05), 0.);
   BOOST_CHECK_CLOSE(scalar_to_fermion_fermion(100., 40., 30., 0.1, 0.05),
                     0.009860002883391398, 1.0e-10);
}
