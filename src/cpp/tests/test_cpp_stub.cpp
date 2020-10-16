//Tests for cpp_stub

#include<cpp_stub.h>
#include<sstream>
#include<fstream>

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

int testSayHello(std::ofstream &results){
    /*!
     * Test message printed to stdout in sayHello function
     *
     * :param std::ofstream &results: The output file
     */

    std::stringbuf buffer;
    cout_redirect rd(&buffer);

    std::string message = "World!";
    cppStub::sayHello(message);

    std::string result = buffer.str();
    std::string answer = "\nHello World!\n";

    if (result.compare(answer) != 0){
        std::cout << "result.compare( answer ) " << result.compare( answer ) << "\n";
        results << "testFoo & False\n";
        return 1;
    }

    results << "testFoo & True\n";
    return 0;
}

int main(){
    /*!
    The main loop which runs the tests defined in the
    accompanying functions. Each function should output
    the function name followed by & followed by True or False
    if the test passes or fails respectively.
    */

    //Open the results file
    std::ofstream results;
    results.open("results.tex");

    //Run the tests
    testSayHello(results);

    //Close the results file
    results.close();

    return 0;
}
