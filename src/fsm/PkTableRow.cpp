/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/PkTableRow.h"

PkTableRow::PkTableRow(const int maxInput)
	: io(maxInput), i2p(maxInput)
{

}

PkTableRow::PkTableRow(const IOMap & io, const I2PMap & i2p)
	: io(io), i2p(i2p)
{

}

IOMap PkTableRow::getIOMap() const
{
	return io;
}

I2PMap PkTableRow::getI2PMap() const
{
	return i2p;
}

int PkTableRow::get(const int x) const
{
	return i2p.at(x);
}

bool operator==(PkTableRow const & row, S2CMap const & s2c)
{
	for (unsigned int i = 0; i < row.i2p.size(); ++ i)
	{
		if (s2c.at(row.i2p.at(i)) != s2c.at(row.get(i)))
		{
			return false;
		}
	}
	return true;
}

std::ostream & operator<<(std::ostream & out, const PkTableRow & pkTableRow)
{
	for (unsigned int i = 0; i < pkTableRow.i2p.size(); ++ i)
	{
		out << " & " << pkTableRow.get(i);
	}
	out << "\\\\\\hline" << std::endl;
	return out;
}
