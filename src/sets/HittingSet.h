/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_SETS_HITTINGSET_H_
#define FSM_SETS_HITTINGSET_H_

#include <unordered_set>
#include <vector>

class HittingSet
{
private:
	/**
	 * The collection of sets ("set system") for which
     * the minimal hitting set should be calculated
	 */
	std::vector<std::unordered_set<int>> s;

	/** the current candidate for the minimal hitting set problem */
	std::unordered_set<int> h;
public:
   /**
	* Create an object for solving the minimal hitting set problem.
	* @param s The set system, i.e. the
    *          collection of sets for which the minimal
    *          hitting set shall be calculated.
	*/
	HittingSet(const std::vector<std::unordered_set<int>>& s);

	/**
	 * Calculate the smallest hitting set for the set system
     * specified when instantiating the object.
	 * @return The smallest set into the hitting set
     *
     * @note this algorithm has worst case complexity
     * of O(2^(#(union s)))
	 */
	std::unordered_set<int> calcMinCardHittingSet() const;
};
#endif //FSM_SETS_HITTINGSET_H_
