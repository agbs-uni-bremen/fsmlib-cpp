/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/OFSMTableRow.h"

OFSMTableRow::OFSMTableRow(const int maxInput, const int maxOutput)
	: maxInput(maxInput), maxOutput(maxOutput)
{
	matrix.insert(matrix.end(), maxInput + 1, std::vector<int>());
	for (int i = 0; i <= maxInput; ++ i)
	{
		matrix.at(i).insert(matrix.at(i).end(), maxOutput + 1, -1);
	}
}

void OFSMTableRow::set(const int i, const int j, const int postState)
{
	matrix [i] [j] = postState;//insertion
}

int OFSMTableRow::get(const int i, const int j) const
{
	return matrix.at(i).at(j);
}

bool OFSMTableRow::ioEquals(const std::shared_ptr<OFSMTableRow> r) const
{
	for (int i = 0; i <= maxInput; ++ i)
	{
		for (int j = 0; j <= maxOutput; ++ j)
		{
			if ((get(i, j) >= 0 && r->get(i, j) < 0) || (get(i, j) < 0 && r->get(i, j) >= 0))
			{
				return false;
			}
		}
	}
	return true;
}

bool OFSMTableRow::classEquals(const S2CMap & s2c, const std::shared_ptr<OFSMTableRow> r)
{
	if (!ioEquals(r))
	{
		return false;
	}

	for (int i = 0; i <= maxInput; ++ i)
	{
		for (int j = 0; j <= maxOutput; ++ j)
		{
			if (matrix.at(i).at(j) < 0)
			{
				continue;
			}
			if (s2c.at(matrix.at(i).at(j)) != s2c.at(r->get(i, j)))
			{
				return false;
			}
		}
	}
	return true;
}
