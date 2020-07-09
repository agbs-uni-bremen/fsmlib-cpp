/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "trees/TestSuite.h"

#include <iostream>
#include <fstream>

using namespace std;

TestSuite::TestSuite()
{

}

bool TestSuite::isEquivalentTo(TestSuite& theOtherTs,
                               bool writeOutput)
{
	bool pass = true;

	if (size() != theOtherTs.size())
	{
		if ( writeOutput )
            cout << "Test suites have different sizes" << endl;
		return false;
	}

	for (unsigned int i = 0; i < size(); ++ i)
	{
		if ( writeOutput ) cout << "Test Case No. " << i << ": ";
		if (!((*this) [i] == theOtherTs [i]))
		{
			if (writeOutput) {
				cout << "FAIL" << endl << "Discrepancy found." <<
					endl << "Mine = " << (*this)[i] << endl <<
					"Other = " << theOtherTs[i] << endl;
			}
			pass = false;
		}
		else
		{
			if ( writeOutput ) cout << "PASS" << endl;
		}
	}
	return pass;
}

bool TestSuite::isReductionOf(TestSuite& theOtherTs,
                              bool writeOutput)
{
	bool pass = true;

	if (size() != theOtherTs.size())
	{
		cout << "Test suites have different sizes" << endl;
		return false;
	}

	for (unsigned int i = 0; i < size(); ++ i)
	{
		if (!theOtherTs [i].contains((*this) [i]))
		{
			if (writeOutput) {
				cout << "Discrepancy found." << endl << "Mine = " << (*this)[i] << endl << "Other = " << theOtherTs[i] << endl;
			}
			pass = false;
		}
	}
	return pass;
}

ostream & operator<<(ostream & out, const TestSuite & testSuite)
{
	for (OutputTree ot : testSuite)
	{
		out << ot;
	}
	return out;
}

void TestSuite::save(const std::string &name) {
    
    ofstream out(name);
    
    out << *this;
    
    out.close();
    
}

size_t TestSuite::totalLength() const
{
    size_t length = 0;
    for(unsigned i = 0; i < this->size(); ++i)
    {
        OutputTree o = this->at(i);
        length += o.getInputTrace().get().size();
    }
    return length;
}
