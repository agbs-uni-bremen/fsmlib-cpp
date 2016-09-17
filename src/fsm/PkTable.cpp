/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/PkTable.h"
#include "fsm/FsmNode.h"
#include "fsm/Dfsm.h"
#include "fsm/FsmLabel.h"
#include "fsm/FsmTransition.h"

PkTable::PkTable(const int numStates, const int maxInput, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: s2c(numStates), maxInput(maxInput), presentationLayer(presentationLayer)
{
	rows.insert(rows.end(), numStates, nullptr);
}

PkTable::PkTable(const int numStates, const int maxInput, const std::vector<std::shared_ptr<PkTableRow>> rows, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: rows(rows), s2c(numStates), maxInput(maxInput), presentationLayer(presentationLayer)
{

}

void PkTable::setRow(const int s, const std::shared_ptr<PkTableRow> row)
{
	rows [s] = row;//insertion
}

void PkTable::setClass(const int n, const int c)
{
	s2c [n] = c;//insertion
}

int PkTable::getClass(const int n) const
{
	return s2c.at(n);
}

int PkTable::maxClassId() const
{
	int id = 0;
	for (unsigned int i = 0; i < s2c.size(); ++ i)
	{
		if (s2c.at(i) > id)
		{
			id = s2c.at(i);
		}
	}
	return id;
}

std::shared_ptr<PkTable> PkTable::getPkPlusOneTable() const
{
	std::shared_ptr<PkTable> pkp1 = std::make_shared<PkTable>(rows.size(), maxInput, rows, presentationLayer);

	int thisClass = 0;
	int thisNewClassId = maxClassId() + 1;
	std::shared_ptr<PkTableRow> refRow;
	std::shared_ptr<PkTableRow> newClassRefRow;
	bool haveNewClasses = false;

	do
	{
		refRow = nullptr;
		for (unsigned int i = 0; i < rows.size(); ++ i)
		{
			if (s2c.at(i) != thisClass)
			{
				continue;
			}

			if (pkp1->getClass(i) >= 0)
			{
				continue;
			}

			if (refRow == nullptr)
			{
				refRow = rows.at(i);
				pkp1->setClass(i, thisClass);
				continue;
			}

			if (*rows.at (i) == s2c)
			{
				pkp1->setClass(i, thisClass);
			}
			else
			{
				haveNewClasses = true;
				newClassRefRow = rows.at(i);
				pkp1->setClass(i, thisNewClassId);

				for (unsigned int j = i + 1; j < rows.size(); ++ j)
				{
					if (s2c.at(j) == thisClass && *rows.at(j) == s2c)
					{
						pkp1->setClass(j, thisNewClassId);
					}
				}

				newClassRefRow = nullptr;
				++ thisNewClassId;
			}
		}
		++ thisClass;
	} while (refRow != nullptr);

	return haveNewClasses ? pkp1 : nullptr;
}

Dfsm PkTable::toFsm(std::string name)
{
	std::string minFsmName = name + "_MIN";
	std::vector<std::shared_ptr<FsmNode>> nodeLst;

	/*Create the FSM states, one for each class*/
	for (int i = 0; i <= maxClassId(); ++ i)
	{
		std::shared_ptr<FsmNode> newNode = std::make_shared<FsmNode>(i, minFsmName + "\n" + getMembers(i), presentationLayer);
		nodeLst.push_back(newNode);
	}

	/*For each FSM state, add outgoing transitions*/
	for (std::shared_ptr<FsmNode> srcNode : nodeLst)
	{
		int classId = srcNode->getId();
		std::shared_ptr<PkTableRow> row = nullptr;

		for (unsigned int i = 0; i < rows.size() && row == nullptr; ++ i)
		{
			if (classId == s2c.at(i))
			{
				row = rows.at(i);
			}
		}

		for (int x = 0; x <= maxInput; ++ x)
		{
			int y = row->getIOMap().at(x);
			int cAux = row->getI2PMap().at(x);
			int cTarget = s2c.at(cAux);

			std::shared_ptr<FsmNode> tgtNode = nullptr;
			for (std::shared_ptr<FsmNode> node : nodeLst)
			{
				if (node->getId() == cTarget)
				{
					tgtNode = node;
					break;
				}
			}
			FsmLabel lbl = FsmLabel(x, y, presentationLayer);
			srcNode->addTransition(FsmTransition(srcNode, tgtNode, lbl));
		}
	}

	return Dfsm(minFsmName, maxInput, 0, nodeLst, presentationLayer);
}

std::string PkTable::getMembers(const int c) const
{
	std::string memSet = "{";
	bool first = true;
	for (unsigned int i = 0; i < rows.size(); ++ i)
	{
		if (s2c.at(i) != c)
		{
			continue;
		}

		if (!first)
		{
			memSet += ",";
		}
		first = false;
		memSet += std::to_string(i);
	}
	memSet += "}";
	return memSet;
}
