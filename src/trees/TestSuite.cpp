/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "trees/TestSuite.h"

TestSuite::TestSuite()
{

}

bool TestSuite::isEquivalentTo(TestSuite & theOtherTs)
{
	bool pass = true;

	if (size() != theOtherTs.size())
	{
		std::cout << "Test suites have different sizes" << std::endl;
		return false;
	}

	for (unsigned int i = 0; i < size(); ++ i)
	{
		std::cout << "Test Case No. " << i << ": ";
		if (!((*this) [i] == theOtherTs [i]))
		{
			std::cout << "FAIL" << std::endl << "Discrepancy found." << std::endl << "Mine = " << (*this) [i] << std::endl << "Other = " << theOtherTs [i] << std::endl;
			pass = false;
		}
		else
		{
			std::cout << "PASS" << std::endl;
		}
	}
	return pass;
}

bool TestSuite::isReductionOf(TestSuite & theOtherTs)
{
	bool pass = true;

	if (size() != theOtherTs.size())
	{
		std::cout << "Test suites have different sizes" << std::endl;
		return false;
	}

	for (unsigned int i = 0; i < size(); ++ i)
	{
		if (!theOtherTs [i].contains((*this) [i]))
		{
			std::cout << "Discrepancy found." << std::endl << "Mine = " << (*this) [i] << std::endl << "Other = " << theOtherTs [i] << std::endl;
			pass = false;
		}
	}
	return pass;
}

std::ostream & operator<<(std::ostream & out, const TestSuite & testSuite)
{
	for (OutputTree ot : testSuite)
	{
		out << std::endl << ot;
	}
	return out;
}
