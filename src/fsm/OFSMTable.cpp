/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#include "fsm/OFSMTable.h"
#include "fsm/OFSMTableRow.h"
#include "fsm/Fsm.h"
#include "fsm/FsmNode.h"
#include "fsm/FsmTransition.h"

std::shared_ptr<OFSMTable> OFSMTable::nextAfterZero()
{
	std::shared_ptr<OFSMTable> next = std::make_shared<OFSMTable>(numStates, maxInput, maxOutput, rows, presentationLayer);
	next->tblId = 1;

	int thisClass = 0;
	S2CMap newS2C = next->getS2C();
	for (int n = 0; n < numStates; ++ n)
	{
		/*If FSM node n is already associated with a class,
		skip this node*/
		if (newS2C.at(n) >= 0)
		{
			continue;
		}

		/*This node needs a new class*/
		newS2C [n] = thisClass;//insertion

		/*Now check all nodes with id > n, whether they
		have to be associated with the same class.
		(all nodes <= n have their class already defined)*/
		for (int m = n + 1; m < numStates; ++ m)
		{
			/*Node m should be associated with the same class as
			n, if and only if its outgoing transitions are labelled
			in the same way as the ones of node n*/
			if (rows.at(n)->ioEquals(rows.at(m)))
			{
				newS2C [m] = thisClass;//insertion
			}
		}

		/*Now we have captured all nodes that should be assigned with
		thisClass. For the next n that is not yet associated with a class,
		we use the next class number*/
		++ thisClass;
	}

	next->setS2C(newS2C);//not sure if abut this way
	return next;
}

OFSMTable::OFSMTable(const std::vector<std::shared_ptr<FsmNode>>& nodes, const int maxInput, const int maxOutput, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: numStates(static_cast<int> (nodes.size())), maxInput(maxInput), maxOutput(maxOutput), tblId(0), s2c(numStates), presentationLayer(presentationLayer)
{
	for (int n = 0; n < numStates; ++ n)
	{
		s2c [n] = 0;//insertion
	}

	for (int i = 0; i < numStates; ++ i)
	{
		rows.push_back(std::make_shared<OFSMTableRow>(maxInput, maxOutput));

		for (FsmTransition tr : nodes.at(i)->getTransitions())
		{
			int x = tr.getLabel().getInput();
			int y = tr.getLabel().getOutput();
			std::shared_ptr<FsmNode> tgtNode = tr.getTarget();
			rows.at(i)->set(x, y, tgtNode->getId());
		}
	}
}

OFSMTable::OFSMTable(const int numStates, const int maxInput, const int maxOutput, const std::vector<std::shared_ptr<OFSMTableRow>>& rows, const std::shared_ptr<FsmPresentationLayer> presentationLayer)
	: numStates(numStates), maxInput(maxInput), maxOutput(maxOutput), tblId(0), s2c(numStates), rows(rows), presentationLayer(presentationLayer)
{

}

int OFSMTable::getId()
{
	return tblId;
}

S2CMap OFSMTable::getS2C()
{
	return s2c;
}

void OFSMTable::setS2C(const S2CMap & ps2c)
{
	s2c = ps2c;
}

int OFSMTable::get(const int id, const int x, const int y)
{
	return rows.at(id)->get(x, y);
}

int OFSMTable::maxClassId() const
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

std::shared_ptr<OFSMTable> OFSMTable::next()
{
	if (tblId == 0)
	{
		return nextAfterZero();
	}

	std::shared_ptr<OFSMTable> next = std::make_shared<OFSMTable>(numStates, maxInput, maxOutput, rows, presentationLayer);
	next->tblId = tblId + 1;

	int thisClass = 0;
	int thisNewClassId = maxClassId() + 1;
	std::shared_ptr<OFSMTableRow> refRow;
	std::shared_ptr<OFSMTableRow> newClassRefRow;
	S2CMap newS2C = next->getS2C();
	bool haveNewClasses = false;

	do
	{
		refRow = nullptr;

		for (int n = 0; n < numStates; ++ n)
		{
			/*We are only interested in nodes n belonging to
			the same class according to s2c*/
			if (s2c.at(n) != thisClass)
			{
				continue;
			}

			/*If n is already associated with a class in the new OFSM table
			nxt, then do not change this*/
			if (newS2C.at(n) >= 0)
			{
				continue;
			}

			/*This applies to the smallest node n that is
			currently associated with thisClass. This node
			always keeps the old class thisClass in the new
			OFSMTable*/
			if (refRow == nullptr)
			{
				refRow = rows.at(n);
				newS2C [n] = thisClass;//insertion
				continue;
			}

			/*This applies to another node that is currently
			also associated with thisClass, and this node
			should remain in thisClass, because its
			post states are all equivalent to the post states
			of the first node associated with this class (this first
			node is represented by refRow).*/
			if (refRow->classEquals(s2c, rows.at(n)))
			{
				newS2C [n] = thisClass;//insertion
				continue;
			}

			/*If we reach this part of the loop, we have found a
			node that was originally associated with the same class
			as the one represented by refRow, but the node has non-equivalent
			post states at one or more outgoing transition.
			This node gets the next unused class id, which is always
			stored in thisNewClassId.*/
			haveNewClasses = true;
			newClassRefRow = rows.at(n);
			newS2C [n] = thisNewClassId;//insertion

			/*Now search for other nodes with id > n that were
			originally in thisClass, but
			should also belong into the new class thisNewClassId*/
			for (int m = n + 1; m < numStates; ++ m)
			{
				if (s2c.at(m) == thisClass && newClassRefRow->classEquals(s2c, rows.at(m)))
				{
					newS2C [m] = thisNewClassId;//insertion
				}
			}

			newClassRefRow = nullptr;
			++ thisNewClassId;
		}

		/*We are done with thisClass, next we refine
		the set of nodes currently belonging to (thisClass+1),
		if any of those exist.*/
		++ thisClass;
	} while (refRow != nullptr);

	next->setS2C(newS2C);
	return haveNewClasses ? next : nullptr;
}

std::string OFSMTable::getMembers(const int c) const
{
	std::string memSet = "{";
	bool first = true;
	for (int i = 0; i < numStates; ++ i)
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

Fsm OFSMTable::toFsm(const std::string & name) const
{
	std::string minFsmName = name;
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
		/*The node id of the new FsmState corresponds to the
		class id we are currently processing.*/
		int classId = srcNode->getId();

		/*Find the first OFSMTableRow where the associated original FsmState
		belongs to class classId. Since other rows associated with
		the same class have equivalent post-states, we only need to
		find one representative row.*/
		std::shared_ptr<OFSMTableRow> row = nullptr;
		for (int i = 0; i < numStates && row == nullptr; ++ i)
		{
			if (classId == s2c.at(i))
			{
				row = rows [i];
			}
		}

		/*Process all outgoing transitions of the original
		FsmNode represented by row.They also become transitions
		of the new srcNode in the minimised FSM*/
		for (int x = 0; x <= maxInput; ++ x)
		{
			for (int y = 0; y <= maxOutput; ++ y)
			{
				int tgtStateId = row->get(x, y);
				if (tgtStateId >= 0)
				{
					/*get the class id of the target node in the original FSM*/
					int tgtClassId = s2c.at(tgtStateId);

					/*Find the new FsmNode in the minimised FSM
					which has tgtClassId as node id*/
					for (std::shared_ptr<FsmNode> tgtNode : nodeLst)
					{
						if (tgtNode->getId() == tgtClassId)
						{
							/*Create the transition with label x/y
							and target node tgtNode*/
							FsmTransition tr = FsmTransition(srcNode, tgtNode, FsmLabel(x, y, presentationLayer));
							srcNode->addTransition(tr);
							break;
						}
					}
				}
			}
		}
	}
	return Fsm(minFsmName, maxInput, maxOutput, nodeLst, presentationLayer);
}

std::ostream & operator<<(std::ostream & out, const OFSMTable & ofsmTable)
{
	/*Tabular environment: first column for the class, second column
	for the state*/
	out << std::endl;
	out << "\\begin{center}" << std::endl;
	out << "\\begin{tabular}{|c||c||";

	/*One additional column for each input/output combination x/y*/
	for (int x = 0; x <= ofsmTable.maxInput; ++ x)
	{
		for (int y = 0; y <= ofsmTable.maxOutput; ++ y)
		{
			out << "c|";
		}
	}
	out << "}\\hline\\hline" << std::endl;

	/*Table headings: [q] || q || 0/0 | 0/1 | ... | maxInput/maxOutput*/
	out << "{\\bf [q]} & {\\bf q}";
	for (int x = 0; x <= ofsmTable.maxInput; ++ x)
	{
		for (int y = 0; y <= ofsmTable.maxOutput; ++ y)
		{
			out << " & " << x << "/" << y;
		}
	}
	out << "\\\\\\hline\\hline" << std::endl;

	/*Print the table contents*/
	for (int i = 0; i < ofsmTable.numStates; ++ i)
	{
		out << ofsmTable.s2c.at(i) << " & " << i;

		for (int x = 0; x <= ofsmTable.maxInput; ++ x)
		{
			for (int y = 0; y <= ofsmTable.maxOutput; ++ y)
			{
				out << " & " << ofsmTable.rows [i]->get(x, y);
			}
		}
		out << "\\\\\\hline" << std::endl;
	}

	out << "\\hline" << std::endl << "\\end{tabular}" << std::endl << "\\end{center}" << std::endl << std::endl;
	return out;
}