/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/DFSMTableRow.h"

DFSMTableRow::DFSMTableRow(const int q, const int maxInput)
	: state(q), ioSection(maxInput), i2postSection(maxInput)
{

}

IOMap& DFSMTableRow::getioSection()
{
	return ioSection;
}

I2PMap& DFSMTableRow::geti2postSection()
{
	return i2postSection;
}

std::ostream & operator<<(std::ostream & out, const DFSMTableRow & dfsmTableRow)
{
	out << dfsmTableRow.state;

	for (unsigned int i = 0; i < dfsmTableRow.ioSection.size(); ++ i)
	{
		out << " & " << dfsmTableRow.ioSection.at(i);
	}

	for (unsigned int i = 0; i < dfsmTableRow.i2postSection.size(); ++ i)
	{
		out << " & " << dfsmTableRow.i2postSection.at(i);
	}
	out << "\\\\\\hline" << std::endl;
	return out;
}
