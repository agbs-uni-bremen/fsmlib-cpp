/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "sets/HittingSet.h"

HittingSet::HittingSet(const std::vector<std::unordered_set<int>>& s)
	: s(s)
{
    // The initial hitting set candidate h is the union of
    // all sets in s
	for (std::unordered_set<int> z : s)
	{
		h.insert(z.begin(), z.end());
	}
	HsTreeNode::hSmallest = h;
}

std::unordered_set<int> HittingSet::calcMinCardHittingSet() const
{
	HsTreeNode root = HsTreeNode(h, s);
	root.expandNode();
	return HsTreeNode::hSmallest;
}
