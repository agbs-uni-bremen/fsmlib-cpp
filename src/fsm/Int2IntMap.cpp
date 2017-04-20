/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/Int2IntMap.h"

Int2IntMap::Int2IntMap(const int maxInput)
{
	for (int i = 0; i <= maxInput; ++ i)
	{
		this->insert(std::pair<int, int>(i, -1));
	}
}
