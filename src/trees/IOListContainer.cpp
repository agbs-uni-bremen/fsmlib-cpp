/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "trees/IOListContainer.h"

bool IOListContainer::isLastLst(const int maxInput, const std::vector<int>& lst) const
{
	for (auto it = lst.cbegin(); it != lst.cend(); ++ it)
	{
		if (*it < maxInput)
		{
			return false;
		}
	}
	return true;
}

std::vector<int> IOListContainer::nextLst(const int maxInput, const std::vector<int>& lst) const
{
	std::vector<int> nextl;
	if (isLastLst(maxInput, lst))
	{
		return nextl;
	}

	/*We know that at least one input element in lst can still be incremented.
	Its successors must be set to zero. Beware, it is a descending iterator*/
	for (auto it = lst.rbegin(); it != lst.rend(); ++ it)
	{
		if (*it == maxInput)
		{
			nextl.insert(nextl.begin(), 0);
		}
		else
		{
			/*Insert the incremented value at this place into nextl*/
			nextl.insert(nextl.begin(), (*it) + 1);

			/*Now copy the remaining values from lst into nextl*/
			while (++ it != lst.rend())
			{
				nextl.insert(nextl.begin(), *it);
			}
			return nextl;
		}
	}
	return nextl;
}

IOListContainer::IOListContainer(const std::shared_ptr<std::vector<std::vector<int>>> iolLst, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: iolLst(iolLst), presentationLayer(presentationLayer)
{

}

IOListContainer::IOListContainer(const int maxInput, const int minLength, const int maxLenght, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: iolLst(std::make_shared<std::vector<std::vector<int>>>()), presentationLayer(presentationLayer)
{
	for (int len = minLength; len <= maxLenght; ++ len)
	{
		/*Initial list of length len only contains zeroes*/
		std::vector<int> lst;
		for (int j = 0; j < len; ++ j)
		{
			lst.push_back(0);
		}
		iolLst->push_back(lst);

		for (lst = nextLst(maxInput, lst); !lst.empty(); lst = nextLst(maxInput, lst))
		{
			iolLst->push_back(lst);
		}
	}
}

std::shared_ptr<std::vector<std::vector<int>>> IOListContainer::getIOLists() const
{
	return iolLst;
}

void IOListContainer::add(const Trace & trc)
{
	iolLst->push_back(trc.get());
}

int IOListContainer::size() const
{
	return static_cast<int> (iolLst->size());
}

std::ostream & operator<<(std::ostream & out, const IOListContainer & ot)
{
	out << "{ ";

	bool isFirst = true;
	for (std::vector<int>& iLst : *ot.iolLst)
	{
		if (!isFirst)
		{
			out << "," << std::endl << "  ";
		}

		for (unsigned int i = 0; i < iLst.size(); ++ i)
		{
			if (i > 0)
			{
				out << ".";
			}

			if (iLst.at(i) == -1)
			{
				out << "eps";
			}
			else
			{
				out << ot.presentationLayer->getInId(iLst.at(i));
			}
		}
		isFirst = false;
	}
	out << " }";
	return out;
}
