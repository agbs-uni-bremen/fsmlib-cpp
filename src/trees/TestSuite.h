/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_TREES_TESTSUITE_H_
#define FSM_TREES_TESTSUITE_H_

#include <ostream>
#include <vector>

class OutputTree;

class TestSuite : public std::vector<OutputTree>
{
public:
   /**
	* Create a new empty test suite
	*/
	TestSuite();

	/**
	 * Check whether or not this test suite is equivalent to an other one
	 * @param theOtherTs The test suite to compare with this one
     * @param writeOutput if true, the method will write PASS/FAIL
     *        information including discrepancies between expected
     *        and observed I/O-traces to cout.
	 * @return true if they are the same, false otherwise
	 */
	bool isEquivalentTo(TestSuite& theOtherTs,
                        bool writeOutput = false);

	/**
	 * Check whether or not this test suite is a reduction of an other one
	 * @param theOtherTs The test suite to compare with this one
     * @param writeOutput if true, the method will write PASS/FAIL
     *        information including discrepancies between expected
     *        and observed I/O-traces to cout.
	 * @return true if they are the other test suite contain this one, false otherwise
	*/
    bool isReductionOf(TestSuite& theOtherTs,
                       bool writeOutput = false);
    
    

	/**
	Output the TestSuite to a standard output stream
	@param out The standard output stream to use
	@param testSuite The TestSuite to print
	@return The standard output stream used, to allow user to cascade <<
	*/
	friend std::ostream & operator<<(std::ostream & out, const TestSuite & testSuite);
    
    /**
     *   Save test suite to file, using the forma of the << operator
     */
    void save(const std::string &name);

    /**
     * @return The sum of all test case sizes in this test suite
     */
    size_t totalLength() const;
};
#endif //FSM_TREES_TESTSUITE_H_
