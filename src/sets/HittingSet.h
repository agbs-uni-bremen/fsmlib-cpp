/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_SETS_HITTINGSET_H_
#define FSM_SETS_HITTINGSET_H_

#include <iostream>
#include <unordered_set>
#include <vector>

#include "sets/HsTreeNode.h"

class HittingSet
{
private:
	/**
	The union of sets into this hitting set
	*/
	std::vector<std::unordered_set<int>> s;

	//TODO
	std::unordered_set<int> h;
public:
	/**
	Create a new Hitting set
	@param s The union of sets to be inserted into the hitting set
	*/
	HittingSet(const std::vector<std::unordered_set<int>>& s);

	/**
	Calculate the the smallest set
	@return The smallest set into the hitting set
	*/
	std::unordered_set<int> calcMinCardHittingSet() const;
};
#endif //FSM_SETS_HITTINGSET_H_
