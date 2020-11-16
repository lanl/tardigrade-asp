/**
  * \file test_cpp_stub.cpp
  *
  * Tests for cpp_stub
  */

#include<cpp_stub.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_cpp_stub
#include <boost/test/included/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer)
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

BOOST_AUTO_TEST_CASE( testSayHello ){
    /*!
     * Test message printed to stdout in sayHello function
     *
     * :param std::ofstream &results: The output file
     */

    std::stringbuf buffer;
    cout_redirect rd(&buffer);

    std::string message = "World!";
    boost::test_tools::output_test_stream result; 
    {
        cout_redirect guard( result.rdbuf() );
        cppStub::sayHello(message);
    }

    std::string answer = "Hello World!\n";

    BOOST_CHECK( result.is_equal( answer ) );

}
