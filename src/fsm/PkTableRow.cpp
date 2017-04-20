/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/PkTableRow.h"


PkTableRow::PkTableRow(IOMap& io,I2PMap& i2p)
	: io(io), i2p(i2p)
{

}

IOMap& PkTableRow::getIOMap() const
{
	return io;
}

I2PMap& PkTableRow::getI2PMap() const
{
	return i2p;
}

int PkTableRow::get(const int x) const
{
	return i2p.at(x);
}

bool PkTableRow::isEquivalent(const PkTableRow& row, const S2CMap& s2c)
{
	for (unsigned int i = 0; i < i2p.size(); ++i)
	{
		if ( s2c.at(i2p.at(i)) != s2c.at(row.getI2PMap().at(i)) )
        {
			return false;
		}
	}
	return true;
}

std::ostream & operator<<(std::ostream & out, const PkTableRow & pkTableRow)
{
    
    for (unsigned int i = 0; i < pkTableRow.io.size(); ++i)
    {
        out << " & " << pkTableRow.io.at(i);
    }
    
	for (unsigned int i = 0; i < pkTableRow.i2p.size(); ++i)
	{
		out << " & " << pkTableRow.get(i);
	}
	out << "\\\\\\hline" << std::endl;
	return out;
    
}
