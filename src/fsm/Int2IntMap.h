/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_INT2INTMAP_H_
#define FSM_FSM_INT2INTMAP_H_

#include <map>

/*This class is really smaller than in Java because it inherit from std::map*/
class Int2IntMap : public std::map<int, int>
{
public:
	/**
	Create a standard map and initialise every element to -1
	\param maxInput The size of the map
	*/
	Int2IntMap(const int maxInput);
};
#endif //FSM_FSM_INT2INTMAP_H_