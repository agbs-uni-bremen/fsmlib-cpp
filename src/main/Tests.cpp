//#include <iostream>
//#include <fstream>
#include <memory>
//#include <stdlib.h>
//#include <interface/FsmPresentationLayer.h>
#include <fsm/Dfsm.h>
#include <fsm/Fsm.h>
#include <fsm/FsmNode.h>
//#include <fsm/IOTrace.h>
//#include <fsm/FsmPrintVisitor.h>
//#include <fsm/FsmSimVisitor.h>
//#include <fsm/FsmOraVisitor.h>
//#include <trees/IOListContainer.h>
#include <trees/OutputTree.h>
#include <trees/Tree.h>
//#include <trees/TestSuite.h>
//#include "json/json.h"
//
//#include "sets/HittingSet.h"
//#include "sets/HsTreeNode.h"
//#include <algorithm>
//#include <cmath>
//#include "fsm/PkTableRow.h"
//#include "fsm/PkTable.h"
//#include "fsm/DFSMTable.h"
//#include "fsm/DFSMTableRow.h"
//#include "fsm/OFSMTableRow.h"
//#include "fsm/OFSMTable.h"
#include "fsm/FsmTransition.h"
//#include "fsm/FsmLabel.h"

//#include <tuple>

#include "Tests.h"

#include <typeinfo>


using namespace std;
//using namespace Json;


// Selects two nodes of m randomly and changes the transitions of one of those nodes in a way that makes both equivalent.
// It is expected that srand() was called before.
shared_ptr<Fsm> makeStatesEquivalent(const Fsm &m) {
	vector<shared_ptr<FsmNode> > lst;
	for (int n = 0; n <= m.getMaxState(); n++) {
		lst.push_back(make_shared<FsmNode>(n, m.getName(), m.getPresentationLayer()));
	}

	unsigned int first = rand() % m.getNodes().size();
	unsigned int second = rand() % m.getNodes().size();
	if (first == second) second = (second + 1) % m.getNodes().size();

	// Now add transitions that correspond exactly to the transitions in m,
	for (int n = 0; n <= m.getMaxState(); n++) {
		auto theNewFsmNodeSrc = lst[n];
		auto theOldFsmNodeSrc = m.getNodes()[n];
		for (auto tr : theOldFsmNodeSrc->getTransitions()) {
			int tgtId = tr->getTarget()->getId();
			auto newLbl = make_shared<FsmLabel>(*(tr->getLabel()));
			shared_ptr<FsmTransition> newTr =
				make_shared<FsmTransition>(theNewFsmNodeSrc, lst[tgtId], newLbl);
			theNewFsmNodeSrc->addTransition(newTr);
		}
	}
	// replaced transitions of lst[second] by transitions with the same labels and targets as the transitions of 
	// lst[first]. This results in two equivalent nodes iff first != second.
	if (first != second) {
		lst[second]->getTransitions().clear();
		for (const auto tr : lst[first]->getTransitions()) {
			lst[second]->getTransitions().push_back(make_shared<FsmTransition>(lst[second], tr->getTarget(), make_shared<FsmLabel>(*(tr->getLabel()))));
		}
	}
	int mI = 0;
	int mO = 0;
	for (const auto n : lst) {
		for (const auto tr : n->getTransitions()) {
			if (tr->getLabel()->getInput() > mI) mI = tr->getLabel()->getInput();
			if (tr->getLabel()->getOutput() > mO) mO = tr->getLabel()->getOutput();
		}
	}
	return make_shared<Fsm>(m.getName(), mI, mO, lst, m.getPresentationLayer());
}

// Selects some node of m at random and deletes all incoming transitions of this node in the new Fsm.
// The returned Fsm contains unreachable nodes if the selected node is not the initial node.
// It is expected that srand() was called before.
shared_ptr<Fsm> makeStatesUnreachable(const Fsm &m) {
	vector<shared_ptr<FsmNode> > lst;
	for (int n = 0; n <= m.getMaxState(); n++) {
		lst.push_back(make_shared<FsmNode>(n, m.getName(), m.getPresentationLayer()));
	}

	int nodeIdx = rand() % m.getNodes().size();

	// Now add transitions that correspond exactly to the transitions in m, but ignore transitions
	// to node at nodeIdx
	for (int n = 0; n <= m.getMaxState(); n++) {
		auto theNewFsmNodeSrc = lst[n];
		auto theOldFsmNodeSrc = m.getNodes()[n];
		for (auto tr : theOldFsmNodeSrc->getTransitions()) {
			int tgtId = tr->getTarget()->getId();
			if (tgtId == nodeIdx) continue;
			auto newLbl = make_shared<FsmLabel>(*(tr->getLabel()));
			shared_ptr<FsmTransition> newTr =
				make_shared<FsmTransition>(theNewFsmNodeSrc, lst[tgtId], newLbl);
			theNewFsmNodeSrc->addTransition(newTr);
		}
	}
	int mI = 0;
	int mO = 0;
	for (const auto n : lst) {
		for (const auto tr : n->getTransitions()) {
			if (tr->getLabel()->getInput() > mI) mI = tr->getLabel()->getInput();
			if (tr->getLabel()->getOutput() > mO) mO = tr->getLabel()->getOutput();
		}
	}
	return make_shared<Fsm>(m.getName(), mI, mO, lst, m.getPresentationLayer());
}

// Selects some input in each state of m and removes each transition for this input.
// It is expected that srand() was called before.
shared_ptr<Fsm> makeStatesPartial(const shared_ptr<Fsm> m) {
	vector<shared_ptr<FsmNode> > lst;
	for (int n = 0; n <= m->getMaxState(); n++) {
		lst.push_back(make_shared<FsmNode>(n, m->getName(), m->getPresentationLayer()));
	}

	// Now add transitions that correspond exactly to the transitions in
	// m, but ignore transitions with some randomly selected input
	for (int n = 0; n <= m->getMaxState(); n++) {
		auto theNewFsmNodeSrc = lst[n];
		auto theOldFsmNodeSrc = m->getNodes()[n];
		int ignoreInput = rand() % (m->getMaxInput() + 1);
		for (auto tr : theOldFsmNodeSrc->getTransitions()) {
			if (tr->getLabel()->getInput() == ignoreInput) continue;
			int tgtId = tr->getTarget()->getId();
			auto newLbl = make_shared<FsmLabel>(*(tr->getLabel()));
			shared_ptr<FsmTransition> newTr =
				make_shared<FsmTransition>(theNewFsmNodeSrc, lst[tgtId], newLbl);
			theNewFsmNodeSrc->addTransition(newTr);
		}
	}
	int mI = 0;
	int mO = 0;
	for (const auto n : lst) {
		for (const auto tr : n->getTransitions()) {
			if (tr->getLabel()->getInput() > mI) mI = tr->getLabel()->getInput();
			if (tr->getLabel()->getOutput() > mO) mO = tr->getLabel()->getOutput();
		}
	}
	return make_shared<Fsm>(m->getName(), mI, mO, lst, m->getPresentationLayer());
}



void fsmlib_assert(string tc, bool verdict, string comment = "");


void fsmlib_assert(string tc, bool verdict, bool &accPass, string comment = "") {
	accPass = accPass and verdict;
	fsmlib_assert(tc, verdict, comment);
}

/*
	Calculates the set of transitions labels of outgoing transitions from given nodes set.
*/
unordered_set<FsmLabel> calcLblSet(std::unordered_set < shared_ptr<FsmNode>> &nodes) {
	unordered_set<FsmLabel> lblSet;
	for (auto n : nodes) {
		for (auto tr : n->getTransitions()) {
			lblSet.insert(*tr->getLabel());
		}
	}
	return lblSet;
}

template<typename T>
bool containsPair(T &lst, std::pair<std::unordered_set<std::shared_ptr<FsmNode>>, std::unordered_set<std::shared_ptr<FsmNode>>> &pair)
{
	typename T::const_iterator it;
	for (it = lst.begin(); it != lst.end(); ++it)
	{
		if (it->first == pair.first && it->second == pair.second) {
			return true;
		}
	}
	return false;
}



/*
	Calculates and returns the set of target nodes reached from states contained in given node list with transitions labeled with given lbl.
*/
unordered_set<shared_ptr<FsmNode>> calcTargetNodes(unordered_set<shared_ptr<FsmNode>> &nodes, FsmLabel &lbl) {
	unordered_set<shared_ptr<FsmNode>> tgtNds;
	for (auto n : nodes) {
		for (auto tr : n->getTransitions()) {
			if (*tr->getLabel() == lbl) tgtNds.insert(tr->getTarget());
		}
	}
	return tgtNds;
}

/**
	This function returns true iff L(q) = L(u). Otherwise the function returns false.
*/
bool ioEquivalenceCheck(const std::shared_ptr<FsmNode> q, const std::shared_ptr<FsmNode> u) {
	// init wl with (q,u)
	std::deque<std::pair<std::unordered_set<std::shared_ptr<FsmNode>>, std::unordered_set<std::shared_ptr<FsmNode>>>> wl;
	std::unordered_set<std::shared_ptr<FsmNode>> l = std::unordered_set<std::shared_ptr<FsmNode>>{ q };
	std::unordered_set<std::shared_ptr<FsmNode>> r = std::unordered_set<std::shared_ptr<FsmNode>>{ u };
	//wl.push_back(std::pair<std::unordered_set<std::shared_ptr<FsmNode>>, std::unordered_set<std::shared_ptr<FsmNode>>>(first, second));
	wl.push_back({ l,r });
	std::vector<std::pair<std::unordered_set<std::shared_ptr<FsmNode>>, std::unordered_set<std::shared_ptr<FsmNode>>>> processed;

	while (not wl.empty()) {
		auto front = wl.front();
		wl.pop_front();
		//calculate lblSet_l
		unordered_set<FsmLabel> lblSet_l = calcLblSet(front.first);
		//calculate lblSet_r
		unordered_set<FsmLabel> lblSet_r = calcLblSet(front.second);
		if (lblSet_l != lblSet_r) return false;

		// insert front to processed if processed does not contain front already
		if (not containsPair(processed, front)) {
			processed.push_back(front);
		}

		for (auto lbl : lblSet_l) {
			// calc tgtNds_l
			unordered_set<shared_ptr<FsmNode>> tgtNds_l = calcTargetNodes(front.first, lbl);
			// calc tgtNds_r
			unordered_set<shared_ptr<FsmNode>> tgtNds_r = calcTargetNodes(front.second, lbl);
			pair<unordered_set<shared_ptr<FsmNode>>, unordered_set<shared_ptr<FsmNode>>> pair{ tgtNds_l, tgtNds_r };
			if (not containsPair(wl, pair) and not containsPair(processed, pair)) {
				wl.push_back(pair);
			}
		}
	}
	return true;
}



/*
	Checks if given fsm contains a pair of states with the same language. Returns true iff fsm contains states q, q' with
	q != q' and L(q) = L(q').
*/
bool hasEquivalentStates(const Fsm &fsm) {
	for (size_t c1 = 0; c1 < fsm.getNodes().size(); c1++) {
		for (size_t c2 = c1 + 1; c2 < fsm.getNodes().size(); c2++) {
			if (ioEquivalenceCheck(fsm.getNodes().at(c1), fsm.getNodes().at(c2))) {
				return true;
			}
		}
	}
	return false;
}

/*
	Check if fsm1 and fsm2 have the same structure (same labeled transitions between nodes with the same indices).
	This method can be used to check if some Fsms structure was changed by some method call.
	Its faster than checking for isomorphism, because of restrictiveness.
	In other words this function tests if fsm1 and fsm2 have identical nodes lists.
*/
bool checkForEqualStructure(const Fsm &fsm1, const Fsm &fsm2) {
	// fsm1 and fsm2 need to be the same size
	if (fsm1.getNodes().size() != fsm2.getNodes().size()) return false;
	for (size_t i = 0; i < fsm1.getNodes().size(); ++i) {
		// each node should have the same number of transitions
		if (fsm1.getNodes().at(i)->getTransitions().size() != fsm2.getNodes().at(i)->getTransitions().size()) return false;
		for (size_t j = 0; j < fsm1.getNodes().at(i)->getTransitions().size(); ++j) {
			auto fsm1Tr = fsm1.getNodes().at(i)->getTransitions().at(j);
			auto fsm2Tr = fsm2.getNodes().at(i)->getTransitions().at(j);
			// compare fsm1Tr and fsm2Tr
			if (fsm1Tr->getSource()->getId() != fsm2Tr->getSource()->getId()
				|| fsm1Tr->getTarget()->getId() != fsm2Tr->getTarget()->getId()
				|| not (*fsm1Tr->getLabel() == *fsm2Tr->getLabel())) {
				return false;
			}
		}
	}
	return true;
}

/*
	Returns a set of all the states of the given Fsm that are reachable from the initial state of that Fsm.
*/
unordered_set<shared_ptr<FsmNode>> getReachableStates(const Fsm &fsm) {
	unordered_set<shared_ptr<FsmNode>> reached{ fsm.getInitialState() };
	deque<shared_ptr<FsmNode>> wl{ fsm.getInitialState() };
	while (not wl.empty()) {
		shared_ptr<FsmNode> q = wl.front();
		wl.pop_front();
		for (auto tr : q->getTransitions()) {
			// add reached target to wl if this target wasn't reached before
			if (reached.insert(tr->getTarget()).second) {
				wl.push_back(tr->getTarget());
			}
		}
	}
	return reached;
}

/*
	Returns true iff fsm.getNodes() contains any of the given node pointers in nodes.
*/
bool contains(const Fsm &fsm, const unordered_set<shared_ptr<FsmNode>> nodes) {
	for (auto n : nodes) {
		if (fsm.contains(n)) return true;
	}
	return false;
}

/*
	Checks if the given fsm is initial connected.
*/
bool isInitialConnected(const Fsm &fsm) {
	auto reachable = getReachableStates(fsm);
	auto nodes = fsm.getNodes();
	unordered_set < shared_ptr<FsmNode> >nodeSet(nodes.cbegin(), nodes.cend());
	return reachable == nodeSet;
}

// =====================================================================================================================================

// Test functions defined for implemented FSM transformations

/**
 * Test function: Fsm::removeUnreachableNodes()
 */
bool testRemoveUnreachableNodes(Fsm &m1, const string &tcID) {
	bool pass = true;
	// get copy of m1
	Fsm copyOfM1 = Fsm(m1);
	vector<shared_ptr<FsmNode>> unreachableNodes;

	// determine set of unreachable nodes in m1
	auto reachable = getReachableStates(m1);
	unordered_set<shared_ptr<FsmNode>> unreachable;
	for (auto n : m1.getNodes()) {
		if (reachable.count(n) == 0) unreachable.insert(n);
	}	

	// use algorithm to transform m1
	bool b = m1.removeUnreachableNodes(unreachableNodes);

	// first check invariant of m1
	bool invariantViolation = not m1.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Invariant holds for M1 after transformation");
	// stop test execution at this point if invariant of m does not hold anymore
	if (invariantViolation) return pass;

	// check properties of m1	
	fsmlib_assert(tcID, isInitialConnected(m1), pass, "Result of removeUnreachableNodes() is initial connected");
	fsmlib_assert(tcID, not contains(m1, unreachable), pass, "Resulting FSM of removeUnreachableNodes() contains none of the nodes that were unreachable before.");

	// check if L(m1) = L(copyOfM1)
	fsmlib_assert(tcID, ioEquivalenceCheck(m1.getInitialState(), copyOfM1.getInitialState()), pass, "removeUnreachableNodes() does not change language of the FSM");

	unordered_set<shared_ptr<FsmNode>> unreachableNodesSet{ unreachableNodes.cbegin(), unreachableNodes.cend() };
	// check b and unreachableNodes
	fsmlib_assert(tcID, (b and (not unreachable.empty())) || (not b and unreachable.empty()), pass, "removeUnreachableNodes() returns true iff FSM contains some unreachable node");
	fsmlib_assert(tcID, (unreachableNodes.size() == unreachable.size()) and (unreachableNodesSet == unreachable), pass, "unreachableNodes contains each unreachable node that was removed");

	return pass;
}

/**
 * Test function: Fsm::transformToObservableFSM()
 */
bool testTransformToObservableFSM(Fsm &m1, const string &tcID) {
	bool pass = true;
	// get copy of m1
	Fsm copyOfM1 = Fsm(m1);

	// use algorithm to transform m1
	Fsm m2 = m1.transformToObservableFSM();

	// first check invariant of m2
	bool invariantViolation = not m2.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Invariant holds for M2 after transformation");
	// stop test execution at this point if invariant of m2 does not hold anymore
	if (invariantViolation) return pass;

	// check properties of m2
	fsmlib_assert(tcID, m2.isObservable(), pass, "M2 is observable after transformToObservable()");

	// check if L(m1) = L(m2)
	fsmlib_assert(tcID, ioEquivalenceCheck(copyOfM1.getInitialState(), m2.getInitialState()), pass, "transformToObservableFSM() does not change the language");

	// check invariant of m1
	invariantViolation = not m1.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Invariant holds for M1 after transformation");
	// stop test execution at this point if invariant of m1 does not hold anymore
	if (invariantViolation) return pass;

	fsmlib_assert(tcID, checkForEqualStructure(m1, copyOfM1), pass, "M1 was not changed by algorithm");

	return pass;
}

/**
 * Test function: Dfsm::minimise()
 */
bool testMinimise_Dfsm(Dfsm &m1, const string &tcID) {
	bool pass = true;
	// get copy of m1
	const Dfsm copyOfM1 = m1;

	// use algorithm to transform m1
	Dfsm m2 = m1.minimise();

	// first check invariant of m2
	bool invariantViolation = not m2.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Invariant holds for M2 after transformation");
	// stop test execution at this point if invariant of m2 does not hold anymore
	if (invariantViolation) return pass;

	// check properties of m2
	fsmlib_assert(tcID, m2.isDeterministic(), pass, "M2 is deterministic after minimise()");
	fsmlib_assert(tcID, isInitialConnected(m2), pass, "M2 is initial connected after minimise()");
	fsmlib_assert(tcID, not hasEquivalentStates(m2), pass, "M2 has no equivalent states after minimise()");

	// check if L(m1) = L(m2)
	fsmlib_assert(tcID, ioEquivalenceCheck(copyOfM1.getInitialState(), m2.getInitialState()), pass, "minimise() does not change the language");

	// check invariant of m1
	invariantViolation = not m1.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Invariant holds for M1 after transformation");
	// stop test execution at this point if invariant of m1 does not hold anymore
	if (invariantViolation) return pass;
	fsmlib_assert(tcID, isInitialConnected(m1), pass, "M1 is initial connected after minimise()");
	fsmlib_assert(tcID, ioEquivalenceCheck(copyOfM1.getInitialState(), m1.getInitialState()), pass, "Language of M1 was not changed by algorithm");

	return pass;
}

/**
 * Test function: Fsm::minimiseObservableFSM()
 */
bool testMinimiseObservableFSM(Fsm &m1, const string &tcID) {
	bool pass = true;
	// get copy of m1
	Fsm copyOfM1 = Fsm(m1);

	// use algorithm to transform m1
	Fsm m2 = m1.minimiseObservableFSM();

	// first check invariant of m2
	bool invariantViolation = not m2.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Invariant holds for M2 after transformation");
	// stop test execution at this point if invariant of m2 does not hold anymore
	if (invariantViolation) return pass;

	// check properties of m2
	fsmlib_assert(tcID, m2.isObservable(), pass, "M2 is observable after minimiseObservable()");
	fsmlib_assert(tcID, not hasEquivalentStates(m2), pass, "M2 has no equivalent states after minimise()");

	// check if L(m1) = L(m2)
	fsmlib_assert(tcID, ioEquivalenceCheck(copyOfM1.getInitialState(), m2.getInitialState()), pass, "minimiseObservableFSM() does not change the language");

	// check invariant of m1
	invariantViolation = not m1.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Invariant holds for M1 after transformation");
	// stop test execution at this point if invariant of m1 does not hold anymore
	if (invariantViolation) return pass;
	fsmlib_assert(tcID, checkForEqualStructure(m1, copyOfM1), pass, "M1 was not changed by algorithm");
	return pass;
}

/**
 * Test function: Fsm::minimise()
 */
bool testMinimise_Fsm(Fsm &m1, const string &tcID) {
	bool pass = true;
	// get copy of m1
	Fsm copyOfM1 = Fsm(m1);

	// use algorithm to transform m1
	Fsm m2 = m1.minimise();

	// first check invariant of m2
	bool invariantViolation = not m2.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Invariant holds for M2 after transformation");
	// stop test execution at this point if invariant of m2 does not hold anymore
	if (invariantViolation) return pass;

	// check properties of m2	
	fsmlib_assert(tcID, m2.isObservable(), pass, "M2 is observable after minimise()");
	fsmlib_assert(tcID, not hasEquivalentStates(m2), pass, "M2 has no equivalent states after minimise()");
	fsmlib_assert(tcID, isInitialConnected(m2), pass, "M2 is initial connected after minimise()");

	// check if L(m1) = L(m2)
	fsmlib_assert(tcID, ioEquivalenceCheck(copyOfM1.getInitialState(), m2.getInitialState()), pass, "minimise() does not change the language");

	// check invariant of m1
	invariantViolation = not m1.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Invariant holds for M1 after transformation");
	// stop test execution at this point if invariant of m1 does not hold anymore
	if (invariantViolation) return pass;

	// check unexpected sideeffects
	fsmlib_assert(tcID, isInitialConnected(m1), pass, "M1 is initial connected after minimise()");
	fsmlib_assert(tcID, ioEquivalenceCheck(copyOfM1.getInitialState(), m1.getInitialState()), pass, "Language of M1 was not changed by algorithm");

	return pass;
}

template <typename T>
struct FsmTransformationTestCase {
	string id;
	shared_ptr<T> m;
};

template <typename T>
void parseFsmTransformationTSFile(const string &testSuitePath, vector<FsmTransformationTestCase<T>> &testSuite) {
	string fname = testSuitePath;
	//vector<FsmTransformationTestCase> testSuite;
	ifstream inputFile(fname);
	if (inputFile.is_open())
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		string line;
		while (getline(inputFile, line))
		{
			vector<string> lineContent;
			stringstream ss(line);
			for (string elem; getline(ss, elem, ';'); lineContent.push_back(elem));
			FsmTransformationTestCase<T> tc;
			tc.id = lineContent.at(0);
			string m1Path = lineContent.at(1);
			if (m1Path == "../../../resources/TestSuites/FSM-Transformations/FSM002.fsm") {
				shared_ptr<FsmNode> n = make_shared<FsmNode>(0, pl);
				//vector<shared_ptr<FsmNode>> lst{ n };
				T m("M", 0, 0, { n }, pl);
				tc.m = make_shared<T>(m);
				testSuite.push_back(tc);
			}
			else {
				shared_ptr<T> m = make_shared<T>(m1Path, pl, "M");
				tc.m = m;
				testSuite.push_back(tc);
			}
		}
		inputFile.close();

	}
	else
	{
		cout << "Unable to open input file" << endl;
		exit(EXIT_FAILURE);
	}
}

/**
 * Test Suite: Fsm::removeUnreachableNodes()
 */
TestResult removeUnreachableNodes_TS() {
	TestResult result("Fsm::removeUnreachableNodes");
	cout << "============================= Start Test of Fsm::removeUnreachableNodes =============================" << endl;
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(59288);
	for (int i = 0; i < 2500; ++i) {
		auto fsm = Fsm::createRandomFsmRepeatable("M", rand() % 4, rand() % 4 + 1, rand() % 10, pl);
		testRemoveUnreachableNodes(*fsm, "TC-Rand-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}
	srand(71);
	for (int i = 0; i < 2500; ++i) {
		auto fsm = makeStatesUnreachable(*Fsm::createRandomFsmRepeatable("M", rand() % 4, rand() % 4 + 1, rand() % 10, pl));
		testRemoveUnreachableNodes(*fsm, "TC-Rand-(MSU)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(MSU)-" + to_string(i));
	}
	srand(1514225);
	for (int i = 0; i < 5000; ++i) {
		Dfsm dfsm("M", rand() % 15 + 1, rand() % 4, rand() % 4 + 1, pl, true);
		testRemoveUnreachableNodes(dfsm, "TC-Rand-(Dfsm)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(Dfsm)-" + to_string(i));
	}

	cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
	shared_ptr<Dfsm> csm = make_shared<Dfsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
	testRemoveUnreachableNodes(*csm, "TC-CSM-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");
	cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
	shared_ptr<Dfsm> fsb = make_shared<Dfsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
	testRemoveUnreachableNodes(*fsb, "TC-FSBC-0") ? ++result.pass : result.fails.push_back("TC-FSBC-0");
	cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
	shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
	testRemoveUnreachableNodes(*gdc, "TC-GDC-0") ? ++result.pass : result.fails.push_back("TC-GDC-0");

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	vector<FsmTransformationTestCase<Fsm>> testSuite;
	parseFsmTransformationTSFile<Fsm>("../../../resources/TestSuites/FSM-Transformations/Fsm_removeUnreachableNodes.testsuite", testSuite);
	for (auto tc : testSuite) {
		testRemoveUnreachableNodes(*tc.m, "TC-Part-" + tc.id) ? ++result.pass : result.fails.push_back("TC-Part-" + tc.id);
	}
	return result;
}

/**
 * Test Suite: Fsm::transformToObservableFSM()
 */
TestResult transformToObservableFSM_TS() {
	TestResult result("Fsm::transformToObservableFSM");
	cout << "============================= Start Test of Fsm::transformToObservableFSM =============================" << endl;
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(84618);
	for (int i = 0; i < 9000; ++i) {
		auto fsm = Fsm::createRandomFsmRepeatable("M1", rand() % 4, rand() % 4 + 1, rand() % 10, pl);
		testTransformToObservableFSM(*fsm, "TC-Rand-" + to_string(i)) 
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}
	srand(119);
	for (int i = 0; i < 1000; ++i) {
		Dfsm dfsm("M", rand() % 15 + 1, rand() % 4, rand() % 4 + 1, pl, true);
		testTransformToObservableFSM(dfsm, "TC-Rand-(Dfsm)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(Dfsm)-" + to_string(i));
	}


	cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
	shared_ptr<Dfsm> csm = make_shared<Dfsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
	testTransformToObservableFSM(*csm, "TC-CSM-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");
	cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
	shared_ptr<Dfsm> fsb = make_shared<Dfsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
	testTransformToObservableFSM(*fsb, "TC-FSBC-0") ? ++result.pass : result.fails.push_back("TC-FSBC-0");
	cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
	shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
	testTransformToObservableFSM(*gdc, "TC-GDC-0") ? ++result.pass : result.fails.push_back("TC-GDC-0");

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	vector<FsmTransformationTestCase<Fsm>> testSuite;
	parseFsmTransformationTSFile<Fsm>("../../../resources/TestSuites/FSM-Transformations/Fsm_transformToObservable.testsuite", testSuite);
	for (auto tc : testSuite) {
		testTransformToObservableFSM(*tc.m, "TC-Part-" + tc.id) ? ++result.pass : result.fails.push_back("TC-Part-" + tc.id);
	}
	return result;
}

/**
 * Test Suite: Dfsm::minimise()
 */
TestResult minimise_Dfsm_TS() {
	TestResult result("Dfsm::minimise");
	cout << "============================= Start Test of Dfsm::minimise =============================" << endl;
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(56958);
	for (int i = 0; i < 10000; ++i) {
		Dfsm dfsm("M", rand() % 15 + 1, rand() % 4, rand() % 4 + 1, pl, true);
		testMinimise_Dfsm(dfsm, "TC-Rand-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}

	cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
	shared_ptr<Dfsm> csm = make_shared<Dfsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
	testMinimise_Dfsm(*csm, "TC-CSM-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");
	cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
	shared_ptr<Dfsm> fsb = make_shared<Dfsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
	testMinimise_Dfsm(*fsb, "TC-FSBC-0") ? ++result.pass : result.fails.push_back("TC-FSBC-0");
	cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
	shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
	testMinimise_Dfsm(*gdc, "TC-GDC-0") ? ++result.pass : result.fails.push_back("TC-GDC-0");

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	vector<FsmTransformationTestCase<Dfsm>> testSuite;
	parseFsmTransformationTSFile<Dfsm>("../../../resources/TestSuites/FSM-Transformations/Dfsm_minimise.testsuite", testSuite);
	for (auto tc : testSuite) {
		testMinimise_Dfsm(*tc.m, "TC-Part-" + tc.id) ? ++result.pass : result.fails.push_back("TC-Part-" + tc.id);
	}
	return result;
}

/**
 * Test Suite: Fsm::minimiseObservableFSM()
 */
TestResult minimiseObservableFSM_TS() {
	TestResult result("Fsm::minimiseObservableFSM");
	cout << "============================= Start Test of Fsm::minimiseObservableFSM =============================" << endl;
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(37580);
	for (int i = 0; i < 2500; ++i) {
		auto fsm = Fsm::createRandomFsmRepeatable("M1", rand() % 4, rand() % 4 + 1, rand() % 10, pl);
		Fsm ofsm = fsm->transformToObservableFSM();
		testMinimiseObservableFSM(ofsm, "TC-Rand-" + to_string(i)) ? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i)) ;
	}
	srand(828);
	for (int i = 0; i < 2500; ++i) {
		auto fsm = makeStatesEquivalent(*Fsm::createRandomFsmRepeatable("M1", rand() % 4, rand() % 4 + 1, rand() % 10, pl));
		Fsm ofsm = fsm->transformToObservableFSM();
		testMinimiseObservableFSM(ofsm, "TC-Rand-(MSE)-" + to_string(i)) ? ++result.pass : result.fails.push_back("TC-Rand-(MSE)-" + to_string(i));
	}
	srand(973);
	for (int i = 0; i < 5000; ++i) {
		Dfsm dfsm("M", rand() % 15 + 1, rand() % 4, rand() % 4 + 1, pl, true);
		testMinimiseObservableFSM(dfsm, "TC-Rand-(Dfsm)-" + to_string(i)) ? ++result.pass : result.fails.push_back("TC-Rand-(Dfsm)-" + to_string(i));
	}

	cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
	shared_ptr<Dfsm> csm = make_shared<Dfsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
	testMinimiseObservableFSM(*csm, "TC-CSM-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");
	cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
	shared_ptr<Dfsm> fsb = make_shared<Dfsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
	testMinimiseObservableFSM(*fsb, "TC-FSBC-0") ? ++result.pass : result.fails.push_back("TC-FSBC-0");
	cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
	shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
	testMinimiseObservableFSM(*gdc, "TC-GDC-0") ? ++result.pass : result.fails.push_back("TC-GDC-0");

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	vector<FsmTransformationTestCase<Fsm>> testSuite;
	parseFsmTransformationTSFile<Fsm>("../../../resources/TestSuites/FSM-Transformations/Fsm_minimiseObservable.testsuite", testSuite);
	for (auto tc : testSuite) {
		testMinimiseObservableFSM(*tc.m, "TC-Part-" + tc.id) ? ++result.pass : result.fails.push_back("TC-Part-" + tc.id);
	}
	return result;
}

/**
 * Test Suite: Fsm::minimise()
 */
TestResult minimise_Fsm_TS() {
	TestResult result("Fsm::minimise");
	cout << "============================= Start Test of Fsm::minimise =============================" << endl;
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(8368);
	for (int i = 0; i < 3334; ++i) {
		auto fsm = Fsm::createRandomFsmRepeatable("M1", rand() % 4, rand() % 4 + 1, rand() % 10, pl);
		testMinimise_Fsm(*fsm, "TC-Rand-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}
	srand(4438);
	for (int i = 0; i < 3333; ++i) {
		auto fsm = makeStatesUnreachable(*Fsm::createRandomFsmRepeatable("M1", rand() % 4, rand() % 4 + 1, rand() % 10, pl));
		testMinimise_Fsm(*fsm, "TC-Rand-(MSU)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(MSU)-" + to_string(i));
	}
	srand(6492);
	for (int i = 0; i < 3333; ++i) {
		auto fsm = makeStatesEquivalent(*Fsm::createRandomFsmRepeatable("M1", rand() % 4, rand() % 4 + 1, rand() % 10, pl));
		testMinimise_Fsm(*fsm, "TC-Rand-(MSE)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(MSE)-" + to_string(i));
	}

	cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
	shared_ptr<Fsm> csm = make_shared<Fsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
	testMinimise_Fsm(*csm, "TC-CSM-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");
	cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
	shared_ptr<Fsm> fsb = make_shared<Fsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
	testMinimise_Fsm(*fsb, "TC-FSBC-0") ? ++result.pass : result.fails.push_back("TC-FSBC-0");
	cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
	shared_ptr<Fsm> gdc = make_shared<Fsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
	testMinimise_Fsm(*gdc, "TC-GDC-0") ? ++result.pass : result.fails.push_back("TC-GDC-0");

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	vector<FsmTransformationTestCase<Fsm>> testSuite;
	parseFsmTransformationTSFile<Fsm>("../../../resources/TestSuites/FSM-Transformations/Fsm_minimise.testsuite", testSuite);
	for (auto tc : testSuite) {
		testMinimise_Fsm(*tc.m, "TC-Part-" + tc.id) ? ++result.pass : result.fails.push_back("TC-Part-" + tc.id);
	}
	return result;
}


// =====================================================================================================================================

// Test functions defined for Fsm::intersect

typedef std::unordered_set<std::shared_ptr<FsmNode>> reachedStates_t;
typedef std::tuple<reachedStates_t, reachedStates_t, reachedStates_t> reachedStatesTuple_t;

template<typename T>
bool containsTuple(T &lst, reachedStatesTuple_t &tuple)
{
	typename T::const_iterator it;
	for (it = lst.begin(); it != lst.end(); ++it)
	{
		if (std::get<0>(*it) == std::get<0>(tuple) && std::get<1>(*it) == std::get<1>(tuple) && std::get<2>(*it) == std::get<2>(tuple)) {
			return true;
		}
	}
	return false;
}

/**
	This function returns true iff L(intersection) = L(m1) ∩ L(m2). Otherwise the function returns false.
*/
bool languageIntersectionCheck(const Fsm &m1, const Fsm &m2, const Fsm &intersection) {
	// init wl with ({m1.q_0},{m2.q_0},{intersection.q_0})	
	deque<reachedStatesTuple_t> wl;
	reachedStates_t a = reachedStates_t{ m1.getInitialState() };
	reachedStates_t b = reachedStates_t{ m2.getInitialState() };
	reachedStates_t c = reachedStates_t{ intersection.getInitialState() };

	wl.push_back({ a,b,c });
	std::vector<reachedStatesTuple_t> processed;

	while (not wl.empty()) {
		auto front = wl.front();
		wl.pop_front();
		//calculate lblSet_a
		unordered_set<FsmLabel> lblSet_a = calcLblSet(std::get<0>(front));
		//calculate lblSet_b
		unordered_set<FsmLabel> lblSet_b = calcLblSet(std::get<1>(front));
		//calculate lblSet_c
		unordered_set<FsmLabel> lblSet_c = calcLblSet(std::get<2>(front));

		// calculate intersection of lblSet_a and lblSet_b
		unordered_set<FsmLabel> lblSet_i;
		for (auto &lbl : lblSet_a) {
			if (lblSet_b.count(lbl) > 0) lblSet_i.insert(lbl);
		}
		if (lblSet_c != lblSet_i) return false;

		// insert front to processed if processed does not contain front already
		if (not containsTuple(processed, front)) {
			processed.push_back(front);
		}

		for (auto lbl : lblSet_c) {
			// calc tgtNds_a
			reachedStates_t tgtNds_a = calcTargetNodes(std::get<0>(front), lbl);
			// calc tgtNds_b
			reachedStates_t tgtNds_b = calcTargetNodes(std::get<1>(front), lbl);
			// calc tgtNds_c
			reachedStates_t tgtNds_c = calcTargetNodes(std::get<2>(front), lbl);
			reachedStatesTuple_t tuple{ tgtNds_a, tgtNds_b, tgtNds_c };
			if (not containsTuple(wl, tuple) and not containsTuple(processed, tuple)) {
				wl.push_back(tuple);
			}
		}
	}
	return true;
}

// Test function for intersect

/**
 * Test function for Fsm::intersect(const Fsm & f).
 */
bool testIntersection(Fsm &m1, const Fsm &m2, const string &tcID) {
	bool pass = true;
	// get copy of m1 and m2
	const Fsm copyOfM1 = Fsm(m1);

	// use Algorithm to calculate result
	const Fsm intersection = m1.intersect(m2);

	// first check invariant for m1 and intersection   (we don't need to check invariant for m2 because it's const)
	bool invariantViolationOfM1 = not m1.checkInvariant();
	fsmlib_assert(tcID, not invariantViolationOfM1, pass, "Invariant still holds for M1 after calculation.");
	bool invariantViolationOfIntersection = not intersection.checkInvariant();
	fsmlib_assert(tcID, not invariantViolationOfIntersection, pass, "Invariant holds for intersection after calculation.");
	// stop test execution at this point if invariant of m or intersection does not hold anymore
	if (invariantViolationOfM1 || invariantViolationOfIntersection) return pass;

	// check language intersection
	fsmlib_assert(tcID, languageIntersectionCheck(m1, m2, intersection), pass, "Language of the result is intersection of L(M1) and L(M2)");

	// check for forbidden side effects
	fsmlib_assert(tcID, checkForEqualStructure(m1, copyOfM1), pass, "M1 was not changed by algorithm");

	return pass;
}

struct IntersectTestCase {
	string id;
	string m1Path;
	string m2Path;
};

shared_ptr<vector<IntersectTestCase>> parseIntersectTSFile(const string &testSuitePath) {
	string fname = testSuitePath;
	vector<IntersectTestCase> testSuite;
	ifstream inputFile(fname);
	if (inputFile.is_open())
	{
		string line;

		while (getline(inputFile, line))
		{
			vector<string> lineContent;
			stringstream ss(line);

			for (string elem; getline(ss, elem, ';'); lineContent.push_back(elem));

			IntersectTestCase tc;
			tc.id = lineContent.at(0);
			tc.m1Path = lineContent.at(1);
			tc.m2Path = lineContent.at(2);

			testSuite.push_back(tc);
		}
		inputFile.close();

		return make_shared<vector<IntersectTestCase>>(testSuite);
	}
	else
	{
		cout << "Unable to open input file" << endl;
		exit(EXIT_FAILURE);
	}
}

/*
 *	Test Suite of Fsm::intersect
*/
TestResult intersect_TS() {
	TestResult result("Fsm::intersect");
	cout << "============================= Start Test of Fsm::intersect =============================" << endl;
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(950);
	for (int i = 0; i < 2500; ++i) {
		auto m1 = Fsm::createRandomFsmRepeatable("M1", rand() % 4, rand() % 4 + 1, rand() % 8, make_shared<FsmPresentationLayer>());
		const auto m2 = m1->createMutantRepeatable("M2", rand() % 5 + 1, rand() % 5 + 1);
		testIntersection(*m1, *m2, "TC-Rand-(1)-" + to_string(i)) ? ++result.pass : result.fails.push_back("TC-Rand-(1)-" + to_string(i));
	}
	srand(104385);
	for (int i = 0; i < 2500; ++i) {
		auto m1 = Fsm::createRandomFsmRepeatable("M1", rand() % 4, rand() % 4 + 1, rand() % 10, make_shared<FsmPresentationLayer>());
		const auto m2 = Fsm::createRandomFsmRepeatable("M2", rand() % 4, rand() % 4 + 1, rand() % 10, make_shared<FsmPresentationLayer>());
		testIntersection(*m1, *m2, "TC-Rand-(2)-" + to_string(i)) ? ++result.pass : result.fails.push_back("TC-Rand-(2)-" + to_string(i));
	}
	srand(158487);
	for (int i = 0; i < 2500; ++i) {
		auto m1 = Dfsm("M1", rand() % 15 + 1, rand() % 4, rand() % 4 + 1, pl, true);
		const auto m2 = m1.createMutantRepeatable("M2", rand() % 5, rand() % 5);
		testIntersection(m1, *m2, "TC-Rand-(3)-" + to_string(i)) ? ++result.pass : result.fails.push_back("TC-Rand-(3)-" + to_string(i));
	}
	srand(627627);
	for (int i = 0; i < 2500; ++i) {
		auto m1 = Dfsm("M1", rand() % 15 + 1, rand() % 4, rand() % 4 + 1, pl, true);
		const auto m2 = Dfsm("M2", rand() % 15 + 1, rand() % 4, rand() % 4 + 1, pl, true);
		testIntersection(m1, m2, "TC-Rand-(4)-" + to_string(i)) ? ++result.pass : result.fails.push_back("TC-Rand-(4)-" + to_string(i));
	}

	srand(353890);
	cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
	shared_ptr<Dfsm> csm = make_shared<Dfsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
	for (int i = 0; i < 100; ++i) {
		testIntersection(*csm, *csm->createMutantRepeatable("Mutant", rand() % 6 + 1, rand() % 6 + 1), "TC-CSM-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-CSM-" + to_string(i));
	}
	cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
	shared_ptr<Dfsm> fsb = make_shared<Dfsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
	for (int i = 0; i < 100; ++i) {
		testIntersection(*fsb, *fsb->createMutantRepeatable("Mutant", rand() % 6 + 1, rand() % 6 + 1), "TC-FSBC-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-FSBC-" + to_string(i));
	}
	cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
	shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
	for (int i = 0; i < 100; ++i) {
		testIntersection(*gdc, *gdc->createMutantRepeatable("Mutant", rand() % 6 + 1, rand() % 6 + 1), "TC-GDC-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-GDC-" + to_string(i));
	}

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	auto testSuite = parseIntersectTSFile("../../../resources/TestSuites/Intersection/Fsm_intersect.testsuite");
	for (auto tc : *testSuite) {
		shared_ptr<Fsm> m1 = make_shared<Fsm>(tc.m1Path, pl, "M1");
		shared_ptr<Fsm> m2 = make_shared<Fsm>(tc.m2Path, pl, "M2");
		testIntersection(*m1, *m2, "TC-Part-" + tc.id) ? ++result.pass : result.fails.push_back("TC-Part-" + tc.id);
	}
	return result;
}


// =====================================================================================================================================

// Test functions for the calculation of distinguishing traces

/*
	Applies inputTrc to startNode and produces the set of outputtraces. Only complete output traces are contained in this
	set (each output trace has the same length as inputTrc).
*/
//set<vector<int>> calcCompleteOutputTraces(const shared_ptr<FsmNode> startNode, const vector<int> inputTrc) {	
//	set<std::tuple<shared_ptr<FsmNode>, vector<int>>> wl{ {startNode, vector<int>()} };
//	set<std::tuple<shared_ptr<FsmNode>, vector<int>>> wl_next;
//
//	for (int x : inputTrc) {
//		for (auto reachedNode : wl) {
//			for (auto transition : std::get<0>(reachedNode)->getTransitions()) {
//				if (transition->getLabel()->getInput() != x) continue; 
//				vector<int> outputTrc = get<1>(reachedNode);
//				outputTrc.push_back(transition->getLabel()->getOutput());
//				wl_next.insert({ transition->getTarget(), outputTrc });
//			}
//		}
//		wl = wl_next;
//		wl_next = set<std::tuple<shared_ptr<FsmNode>, vector<int>>>();
//	}
//
//	set<vector<int>> outputTrcs;
//	for (auto reachedNode : wl) {
//		outputTrcs.insert(std::get<1>(reachedNode));
//	}
//
//	return outputTrcs;
//}

/*
	Second Version of Algorithm. Applies inputTrc to startNode and produces set of outputtraces. If some FsmNode is reached with a prefix
	of inputTrc in which the next input of inputTrc is undefined, the corresponding output trace will be expanded by an 'NULL'/-1 output
	(not contained in the output alphabet) and algorithm stays in this FsmNode. Then the next input is applied.
*/
set<vector<int>> calcCompleteOutputTraces(const shared_ptr<FsmNode> startNode, const vector<int> inputTrc) {
	set<std::tuple<shared_ptr<FsmNode>, vector<int>>> wl{ {startNode, vector<int>()} };
	set<std::tuple<shared_ptr<FsmNode>, vector<int>>> wl_next;
	const int nullOutput = -1;

	for (int x : inputTrc) {
		for (auto reachedNode : wl) {
			bool defined = false;
			for (auto transition : std::get<0>(reachedNode)->getTransitions()) {
				if (transition->getLabel()->getInput() != x) continue;
				defined = true;
				vector<int> outputTrc = get<1>(reachedNode);
				outputTrc.push_back(transition->getLabel()->getOutput());
				wl_next.insert({ transition->getTarget(), outputTrc });
			}
			if (not defined) {
				vector<int> outputTrc = get<1>(reachedNode);
				outputTrc.push_back(nullOutput);
				wl_next.insert({ std::get<0>(reachedNode), outputTrc });
			}
		}
		wl = wl_next;
		wl_next = set<std::tuple<shared_ptr<FsmNode>, vector<int>>>();
	}

	set<vector<int>> outputTrcs;
	for (auto reachedNode : wl) {
		outputTrcs.insert(std::get<1>(reachedNode));
	}
	return outputTrcs;
}

/*
	Second Version of Algorithm. Applies inputTrc to startNode and produces set of outputtraces. If some FsmNode is reached with a prefix
	of inputTrc in which the next input of inputTrc is undefined, the corresponding output trace will be expanded by an 'NULL' output
	(not contained in the output alphabet) and algorithm stays in this FsmNode. Then the next input is applied.
*/
//set<vector<int>> calcCompleteOutputTraces2(const shared_ptr<FsmNode> startNode, const vector<int> inputTrc, const int maxOutput) {
//	set<std::tuple<shared_ptr<FsmNode>, vector<int>>> wl{ {startNode, vector<int>()} };
//	set<std::tuple<shared_ptr<FsmNode>, vector<int>>> wl_next;
//	const int nullOutput = maxOutput + 1;
//
//	for (int x : inputTrc) {
//		for (auto reachedNode : wl) {
//			bool defined = false;
//			for (auto transition : std::get<0>(reachedNode)->getTransitions()) {
//				if (transition->getLabel()->getInput() != x) continue;
//				defined = true;
//				vector<int> outputTrc = get<1>(reachedNode);
//				outputTrc.push_back(transition->getLabel()->getOutput());
//				wl_next.insert({ transition->getTarget(), outputTrc });
//			}
//			if (not defined) {
//				vector<int> outputTrc = get<1>(reachedNode);
//				outputTrc.push_back(nullOutput);
//				wl_next.insert({ std::get<0>(reachedNode), outputTrc });
//			}
//		}
//		wl = wl_next;
//		wl_next = set<std::tuple<shared_ptr<FsmNode>, vector<int>>>();
//	}
//
//	set<vector<int>> outputTrcs;
//	for (auto reachedNode : wl) {
//		//for (auto i : std::get<1>(reachedNode)) cout << i << ",";
//		//cout << "\n";
//		outputTrcs.insert(std::get<1>(reachedNode));
//	}
//	return outputTrcs;
//}

/*
	Checks if given inTrc is a Distinguishing Trace for q1 and q2.
	Returns true iff inTrc produces some outTrc of the same length
	which is only contained in the language of one of these FsmNodes.
*/
bool isDistTrc(const shared_ptr<FsmNode> q1, const shared_ptr<FsmNode> q2, const vector<int> &inTrc) {
	return calcCompleteOutputTraces(q1, inTrc) != calcCompleteOutputTraces(q2, inTrc);	
}

/*
	Returns true iff w contains a Distinguishing Trace for q1 and q2.
*/
bool containsDistTrcForPair(const shared_ptr<FsmNode> q1, const shared_ptr<FsmNode> q2, const IOListContainer &w) {
	for (auto inTrc : *w.getIOLists()) {
		if (isDistTrc(q1, q2, inTrc)) return true;
	}
	return false;
}

/*
	m has to be minimal and observable.
	Returns true iff w is a Characterisation Set of m.
*/
bool isCharaterisationSet(const Fsm &m, const IOListContainer w) {
	for (size_t q1Idx = 0; q1Idx < m.size(); ++q1Idx) {
		for (size_t q2Idx = q1Idx + 1; q2Idx < m.size(); ++q2Idx) {
			if (not containsDistTrcForPair(m.getNodes().at(q1Idx), m.getNodes().at(q2Idx), w)) return false;
		}
	}
	return true;
}

/*
	Checks if trc is a non empty prefix of some element contained in w.
*/
bool isPrefixOfElement(const vector<int> &trc, const std::shared_ptr<const Tree> w) {
	//w->addToRoot() may be faster
	if (trc.size() == 0) return false;

	for (auto elem : *w->getIOLists().getIOLists()) {
		if (elem.size() < trc.size()) continue;

		for (size_t i = 0; i < trc.size(); ++i) {
			if (trc.at(i) != elem.at(i)) continue;
		}
		return true;
	}
	return false;
}

/*
	m has to be minimal and observable. qi ist expected to be a state of m.
	wi is expected to contain traces over the input alphabet of m. w is expected to be a characterisation set of m.

	Returns true iff wi is a State Identification Set for qi in m with Characterisation Set w.
*/
bool isStateIdentificationSet(const Fsm &m, const shared_ptr<FsmNode> qi, const std::shared_ptr<Tree> wi, const std::shared_ptr<Tree> w) {
	for (auto trc : *wi->getIOLists().getIOLists()) {
		if (not isPrefixOfElement(trc, w)) return false;
	}

	for (auto q : m.getNodes()) {
		if (q == qi) continue;

		if (not containsDistTrcForPair(q, qi, wi->getIOLists())) return false;
	}

	return true;
}

/*
	m has to be minimal and observable. qi ist expected to be a state of m.
	wi is expected to be a State Identification Set of qi in m. w is expected to be a characterisation set of m.

	Returns true iff there is no subset of wi that is a State Identification Set of qi in m.
*/
bool isMinimalStateIdentificationSet(const Fsm &m, const shared_ptr<FsmNode> qi, const std::shared_ptr<Tree> wi, const std::shared_ptr<Tree> w) {
	// border case: accept a state identification set with only one element as minimal
	if (wi->getIOLists().getIOLists()->size() == 1) {
		return true;
	}
	auto pl = make_shared<FsmPresentationLayer>();
	for (size_t i = 0; i < wi->getIOLists().getIOLists()->size(); ++i) {
		IOListContainer iolc(pl);
		for (size_t j = 0; j < wi->getIOLists().getIOLists()->size(); ++j) {
			if (i == j) continue;
			iolc.add({ wi->getIOLists().getIOLists()->at(j), pl });
		}
		shared_ptr<Tree> alternativeWi = make_shared<Tree>(make_shared<TreeNode>(), pl);
		alternativeWi->addToRoot(iolc);
		if (isStateIdentificationSet(m, qi, alternativeWi, w)) return false;
	}
	return true;
}



/**
 * Test function for Dfsm::getCharacterisationSet().
 * Parameter m is expected to be a minimal and complete Dfsm.
 */
bool testGetCharacterisationSet_Dfsm(Dfsm &m, const string &tcID) {
	bool pass = true;
	// get copy of m
	const Dfsm copyOfM = Dfsm(m);

	// use Algorithm to calculate result
	const auto w = m.getCharacterisationSet();

	// first check invariant of m
	bool invariantViolation = not m.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Dfsm class invariant still holds for M after calculation.");
	// stop test execution at this point if invariant of m does not hold anymore
	if (invariantViolation) return pass;

	// check definition of 'Characterisation Set' for w
	fsmlib_assert(tcID, isCharaterisationSet(m, w), pass, "Result is a Characterisation Set for M.");

	fsmlib_assert(tcID, *m.characterisationSet->getIOLists().getIOLists() == *w.getIOLists(), pass, "Result is stored in attribute.");

	// check if structure of m has changed
	fsmlib_assert(tcID, checkForEqualStructure(m, copyOfM), pass, "M was not changed by algorithm");

	return pass;
}


/**
 * Test function for Fsm::getCharacterisationSet().
 * Parameter m is expected to be a minimal and observable Fsm.
 */
bool testGetCharacterisationSet_Fsm(Fsm &m, const string &tcID) {
	bool pass = true;
	// get copy of m
	const Fsm copyOfM = Fsm(m);

	// use Algorithm to calculate result
	const auto w = m.getCharacterisationSet();

	// first check invariant of m
	bool invariantViolation = not m.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Fsm class invariant still holds for M after calculation.");
	// stop test execution at this point if invariant of m does not hold anymore
	if (invariantViolation) return pass;

	// check definition of 'Characterisation Set' for w
	fsmlib_assert(tcID, isCharaterisationSet(m, w), pass, "Result is a Characterisation Set for M.");

	fsmlib_assert(tcID, *m.characterisationSet->getIOLists().getIOLists() == *w.getIOLists(), pass, "Result is stored in attribute.");

	// check if structure of m has changed
	fsmlib_assert(tcID, checkForEqualStructure(m, copyOfM), pass, "M was not changed by algorithm");
	return pass;
}

/**
 * Test function for FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
 *                                           const vector<shared_ptr<PkTable>>& pktblLst,
 *                                           const int maxInput)
 * Parameter m is expected to be a minimal and complete Dfsm.
 */
bool testCalcDistinguishingTrace1(Dfsm &m, const string &tcID) {
	bool pass = true;
	// get copy of m
	const Dfsm copyOfM = Dfsm(m);

	// calculate the needed parameters from m
	m.calcPkTables();
	const auto tables = m.pktblLst;

	// test each pair of different nodes
	for (size_t q1Idx = 0; q1Idx < m.size(); ++q1Idx) {
		for (size_t q2Idx = q1Idx + 1; q2Idx < m.size(); ++q2Idx) {
			if (q1Idx == q2Idx) continue;
			const auto q1 = m.getNodes().at(q1Idx);
			const auto q2 = m.getNodes().at(q2Idx);

			// use Algorithm to calculate result
			InputTrace inTrc = q1->calcDistinguishingTrace(q2, tables, m.getMaxInput());

			// first check invariant of m
			bool invariantViolation = not m.checkInvariant();
			fsmlib_assert(tcID, not invariantViolation, pass, "Dfsm class invariant still holds for M after calculation.");
			// stop test execution at this point if invariant of m does not hold anymore
			if (invariantViolation) return pass;

			// check definition of 'Distinguishing Trace' for inTrc
			fsmlib_assert(tcID, isDistTrc(q1, q2, inTrc.get()), pass, "Calculated Trace is a Distinguishing Trace for q1 and q2.");

			// check if structure of m has changed
			fsmlib_assert(tcID, checkForEqualStructure(m, copyOfM), pass, "M was not changed by algorithm");
		}
	}
	return pass;
}

/**
 * Test function for FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
 *                                           const vector<shared_ptr<OFSMTable>>& ofsmTblLst,
 *                                           const int maxInput,
 *                                           const int maxOutput)
 * Parameter m is expected to be a minimal and observable Fsm.
 */
bool testCalcDistinguishingTrace2(Fsm &m, const string &tcID) {
	bool pass = true;
	// get copy of m
	const Fsm copyOfM = Fsm(m);

	// calculate the needed parameters from m
	m.calcOFSMTables();
	const auto tables = m.ofsmTableLst;	

	// test each pair of different nodes
	for (size_t q1Idx = 0; q1Idx < m.size(); ++q1Idx) {
		for (size_t q2Idx = q1Idx + 1; q2Idx < m.size(); ++q2Idx) {
			if (q1Idx == q2Idx) continue;
			const auto q1 = m.getNodes().at(q1Idx);
			const auto q2 = m.getNodes().at(q2Idx);

			// use Algorithm to calculate result
			InputTrace inTrc = q1->calcDistinguishingTrace(q2, tables, m.getMaxInput(), m.getMaxOutput());

			// first check invariant of m
			bool invariantViolation = not m.checkInvariant();
			fsmlib_assert(tcID, not invariantViolation, pass, "Fsm class invariant still holds for M after calculation.");
			// stop test execution at this point if invariant of m does not hold anymore
			if (invariantViolation) return pass;

			// check definition of 'Distinguishing Trace' for inTrc
			fsmlib_assert(tcID, isDistTrc(q1, q2, inTrc.get()), pass, "Calculated Trace is a Distinguishing Trace for q1 and q2.");			

			// check if structure of m has changed
			fsmlib_assert(tcID, checkForEqualStructure(m, copyOfM), pass, "M was not changed by algorithm");
		}
	}
	return pass;
}

/**
 * Test function for Fsm::calcStateIdentificationSets(). 
 * m is expected to be a minimal Fsm.
 */
bool testCalcStateIdentificationSets(Fsm &m, const string &tcID) {
	bool pass = true;
	// get copy of m
	const Fsm copyOfM = Fsm(m);

	// calculate the needed parameters from m
	m.Fsm::getCharacterisationSet();
	const IOListContainer tracesOfW = m.characterisationSet->getIOLists();
	if (tracesOfW.size() > 5) {
		cout << "Characterisation Set is too big. Stop Test Case." << endl;
		return pass;
	}	

	// use Algorithm to calculate result
	m.calcStateIdentificationSets();
	const auto stateIdSets = m.stateIdentificationSets;

	// first check invariant of m
	bool invariantViolation = not m.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Invariant still holds for M after calculation.");
	// stop test execution at this point if invariant of m does not hold anymore
	if (invariantViolation) return pass;

	// Check Definition of minimal State Identification Set for each element in stateIdSets
	fsmlib_assert(tcID, stateIdSets.size() == m.getNodes().size(), pass, "Number of calculated State Identification Sets matches the number of states of M.");
	for (size_t i = 0; i < stateIdSets.size(); ++i) {
		fsmlib_assert(tcID, isStateIdentificationSet(m, m.getNodes().at(i), stateIdSets.at(i), m.characterisationSet), pass, "M.stateIdentificationSets[i] is a State Identification Set for M.nodes[i].");
		fsmlib_assert(tcID, isMinimalStateIdentificationSet(m, m.getNodes().at(i), stateIdSets.at(i), m.characterisationSet), pass, "M.stateIdentificationSets[i] is minimal.");
	}

	// check if structure of m has changed
	fsmlib_assert(tcID, checkForEqualStructure(m, copyOfM), pass, "M was not changed by algorithm");

	// check if m.characterisationSet has changed
	fsmlib_assert(tcID, *tracesOfW.getIOLists() == *m.characterisationSet->getIOLists().getIOLists(), pass, "characterisation set of M has not changed");

	return pass;
}

/**
 * Test function for Fsm::calcStateIdentificationSetsFast().
 */
bool testCalcStateIdentificationSetsFast(Fsm &m, const string &tcID) {
	bool pass = true;
	// get copy of m
	const Fsm copyOfM = Fsm(m);

	// calculate the needed parameters from m
	m.Fsm::getCharacterisationSet();
	const IOListContainer tracesOfW = m.characterisationSet->getIOLists();


	// use Algorithm to calculate result
	m.calcStateIdentificationSetsFast();
	const auto stateIdSets = m.stateIdentificationSets;

	// first check invariant of m
	bool invariantViolation = not m.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Invariant still holds for M after calculation.");
	// stop test execution at this point if invariant of m does not hold anymore
	if (invariantViolation) return pass;

	// Check Definition of State Identification Set for each element in stateIdSets
	fsmlib_assert(tcID, stateIdSets.size() == m.getNodes().size(), "Number of calculated State Identification Sets matches the number of states of M.");
	for (size_t i = 0; i < stateIdSets.size(); ++i) {
		fsmlib_assert(tcID, isStateIdentificationSet(m, m.getNodes().at(i), stateIdSets.at(i), m.characterisationSet), pass, "M.stateIdentificationSets[i] is a State Identification Set for M.nodes[i].");			
	}

	// check if structure of m has changed
	fsmlib_assert(tcID, checkForEqualStructure(m, copyOfM), pass, "M was not changed by algorithm");

	// check if m.characterisationSet has changed
	fsmlib_assert(tcID, *tracesOfW.getIOLists() == *m.characterisationSet->getIOLists().getIOLists(), pass, "characterisation set of M has not changed");
	return pass;
}

template <typename T>
struct DistinguishingTraceTestCase {
	string id;
	shared_ptr<T> m;
};

template <typename T>
void parseDistinguishingTraceTSFile(const string &testSuitePath, vector<DistinguishingTraceTestCase<T>> &testSuite) {
	string fname = testSuitePath;
	ifstream inputFile(fname);
	if (inputFile.is_open())
	{
		shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
		string line;
		while (getline(inputFile, line))
		{
			vector<string> lineContent;
			stringstream ss(line);
			for (string elem; getline(ss, elem, ';'); lineContent.push_back(elem));
			DistinguishingTraceTestCase<T> tc;
			tc.id = lineContent.at(0);
			string mPath = lineContent.at(1);
			tc.m = make_shared<T>(mPath, pl, "M");
			if (lineContent.at(2) == "true") { 
				tc.m->getCharacterisationSet();
			};
			if (lineContent.at(3) == "true") {
				tc.m->calcStateIdentificationSetsFast();
			}; 
			testSuite.push_back(tc);
		}
		inputFile.close();
	}
	else
	{
		cout << "Unable to open input file" << endl;
		exit(EXIT_FAILURE);
	}
}

/*
 *	Random Test Suite for test of Dfsm::getCharacterisationSet().
 */
TestResult getCharacterisationSet_Dfsm_TS() {
	TestResult result("Dfsm::getCharacterisationSet");
	cout << "============================= Start Test of Dfsm::getCharacterisationSet =============================" << endl;

	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	srand(261250);
	for (int i = 0; i < 10000; ++i) {
		auto m = Dfsm("M", rand() % 15 + 1, rand() % 4, rand() % 4 + 1, pl, true);
		auto minM = m.minimise();
		testGetCharacterisationSet_Dfsm(minM, "TC-Rand-" + to_string(i)) ? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}

	cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
	Dfsm csm = Dfsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
	testGetCharacterisationSet_Dfsm(csm, "TC-CSM-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");

	cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
	Dfsm fsb = Dfsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
	testGetCharacterisationSet_Dfsm(fsb, "TC-FSBC-0") ? ++result.pass : result.fails.push_back("TC-FSBC-0");

	cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
	Dfsm gdc = Dfsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
	testGetCharacterisationSet_Dfsm(gdc, "TC-GDC-0") ? ++result.pass : result.fails.push_back("TC-GDC-0");


	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	vector<DistinguishingTraceTestCase<Dfsm>> testSuite;
	parseDistinguishingTraceTSFile("../../../resources/TestSuites/DistinguishingTraces/Dfsm_getCharacterisationSet.testsuite", testSuite);
	for (auto tc : testSuite) {
		testGetCharacterisationSet_Dfsm(*tc.m, "TC-Part-" + tc.id) ? ++result.pass : result.fails.push_back("TC-Part-" + tc.id);
	}
	return result;
}

/*
 *	Random Test Suite for test of Fsm::getCharacterisationSet().
 */
TestResult getCharacterisationSet_Fsm_TS() { 
	TestResult result("Fsm::getCharacterisationSet");
	cout << "============================= Start Test of Fsm::getCharacterisationSet =============================" << endl;
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

	cout << "------------------------------- Start Random Tests -------------------------------" << endl;	
	srand(64162);
	for (int i = 0; i < 5000; ++i) {
		auto fsm = Fsm::createRandomFsmRepeatable("M1", rand() % 4, rand() % 4 + 1, rand() % 6, pl);
		auto minFsm = fsm->minimise();
		string id = "TC-Rand-" + to_string(i);
		testGetCharacterisationSet_Fsm(minFsm, id) ? ++result.pass : result.fails.push_back(id);
	}
	srand(7197);
	for (int i = 0; i < 5000; ++i) {
		auto fsm = makeStatesPartial(Fsm::createRandomFsmRepeatable("M1", rand() % 4 + 1, rand() % 4 + 1, rand() % 6, pl));
		auto minFsm = fsm->minimise();
		string id = "TC-Rand-(MSP)-" + to_string(i);
		testGetCharacterisationSet_Fsm(minFsm, id) ? ++result.pass : result.fails.push_back(id);
	}

	cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
	Fsm csm = Fsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
	testGetCharacterisationSet_Fsm(csm, "TC-CSM-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");

	cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
	Fsm fsb = Fsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
	testGetCharacterisationSet_Fsm(fsb, "TC-FSBC-0") ? ++result.pass : result.fails.push_back("TC-FSBC-0");

	cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
	Fsm gdc = Fsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
	testGetCharacterisationSet_Fsm(gdc, "TC-GDC-0") ? ++result.pass : result.fails.push_back("TC-GDC-0");

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	vector<DistinguishingTraceTestCase<Fsm>> testSuite;
	parseDistinguishingTraceTSFile("../../../resources/TestSuites/DistinguishingTraces/Fsm_getCharacterisationSet.testsuite", testSuite);
	for (auto tc : testSuite) {
		testGetCharacterisationSet_Fsm(*tc.m, "TC-Part-" + tc.id) ? ++result.pass : result.fails.push_back("TC-Part-" + tc.id);
	}
	return result;
}

/*
 *	Random Test Suite for test of FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
 *                                           const vector<shared_ptr<PkTable>>& pktblLst,
 *                                           const int maxInput)
 */
TestResult calcDistinguishingTrace_PkTables_TS() {
	TestResult result("FsmNode::calcDistinguishingTracet(pkTables)");
	cout << "============================= Start Test of FsmNode::calcDistinguishingTrace(pkTables) =============================" << endl;	
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();

	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(329132);
	for (int i = 0; i < 2000; ++i) {
		auto m = Dfsm("M", rand() % 10 + 1, rand() % 4, rand() % 4 + 1, pl, true);
		auto minM = m.minimise();
		testCalcDistinguishingTrace1(minM, "TC-Rand-" + to_string(i)) ? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}

	cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
	Dfsm csm = Dfsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
	testCalcDistinguishingTrace1(csm, "TC-CSM-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");

	cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
	Dfsm fsb = Dfsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
	testCalcDistinguishingTrace1(fsb, "TC-FSBC-0") ? ++result.pass : result.fails.push_back("TC-FSBC-0");

	cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
	Dfsm gdc = Dfsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
	testCalcDistinguishingTrace1(gdc, "TC-GDC-0") ? ++result.pass : result.fails.push_back("TC-GDC-0");

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	vector<DistinguishingTraceTestCase<Dfsm>> testSuite;
	parseDistinguishingTraceTSFile("../../../resources/TestSuites/DistinguishingTraces/FsmNode_calcDistinguishingTrace_pk.testsuite", testSuite);
	for (auto tc : testSuite) {
		testCalcDistinguishingTrace1(*tc.m, "TC-Part-" + tc.id) ? ++result.pass : result.fails.push_back("TC-Part-" + tc.id);
	}
	return result;
}

/*
 *	Random Test Suite for test of FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
 *                                           const vector<shared_ptr<OFSMTable>>& ofsmTblLst,
 *                                           const int maxInput,
 *                                           const int maxOutput)
 */
TestResult calcDistinguishingTrace_OFSMTables_TS() {
	TestResult result("FsmNode::calcDistinguishingTrace(ofsmTables)");
	cout << "============================= Start Test of FsmNode::calcDistinguishingTrace(ofsmTables) =============================" << endl;

	cout << "------------------------------- Start Random Tests -------------------------------" << endl;	
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	srand(36185);
	for (int i = 0; i < 1000; ++i) {
		auto fsm = Fsm::createRandomFsmRepeatable("M1", rand() % 4, rand() % 4 + 1, rand() % 6, pl);
		auto minFsm = fsm->minimise();
		testCalcDistinguishingTrace2(minFsm, "TC-Rand-" + to_string(i)) ? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}
	srand(6075);
	for (int i = 0; i < 1000; ++i) {
		auto fsm = makeStatesPartial(Fsm::createRandomFsmRepeatable("M1", rand() % 4 + 1, rand() % 4 + 1, rand() % 6, pl));
		auto minFsm = fsm->minimise();
		testCalcDistinguishingTrace2(minFsm, "TC-Rand-(MSP)-" + to_string(i)) ? ++result.pass : result.fails.push_back("TC-Rand-(MSP)-" + to_string(i));
	}

	cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
	Fsm csm = Fsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
	testCalcDistinguishingTrace2(csm, "TC-CSM-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");

	cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
	Fsm fsb = Fsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
	testCalcDistinguishingTrace2(fsb, "TC-FSBC-0") ? ++result.pass : result.fails.push_back("TC-FSBC-0");

	cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
	Fsm gdc = Fsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
	testCalcDistinguishingTrace2(gdc, "TC-GDC-0") ? ++result.pass : result.fails.push_back("TC-GDC-0");

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	vector<DistinguishingTraceTestCase<Fsm>> testSuite;
	parseDistinguishingTraceTSFile("../../../resources/TestSuites/DistinguishingTraces/FsmNode_calcDistinguishingTrace_ofsm.testsuite", testSuite);
	for (auto tc : testSuite) {
		testCalcDistinguishingTrace2(*tc.m, "TC-Part-" + tc.id) ? ++result.pass : result.fails.push_back("TC-Part-" + tc.id);
	}
	return result;
}

/*
 *	Random Test Suite for test of Fsm::calcStateIdentificationSets().
 */
TestResult calcStateIdentificationSets_TS() {
	TestResult result("Fsm::calcStateIdentificationSets");
	cout << "============================= Start Test of Fsm::calcStateIdentificationSets() =============================" << endl;

	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	srand(3447);
	for (int i = 0; i < 2500; ++i) {
		auto fsm = Fsm::createRandomFsmRepeatable("M", rand() % 4, rand() % 4 + 1, rand() % 6, make_shared<FsmPresentationLayer>());
		auto minFsm = fsm->minimise();
		cout << minFsm.size() << endl;
		if (minFsm.size() > 50) {
			cout << "M is too big. Stop Test Case." << endl;
			continue;
		}
		testCalcStateIdentificationSets(minFsm, "TC-Rand-" + to_string(i)) 
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}
	srand(76412);
	for (int i = 0; i < 2500; ++i) {
		auto fsm = makeStatesPartial(Fsm::createRandomFsmRepeatable("M", rand() % 4 + 1, rand() % 4 + 1, rand() % 6, make_shared<FsmPresentationLayer>()));
		auto minFsm = fsm->minimise();
		if (minFsm.size() > 50) {
			cout << "M is too big. Stop Test Case." << endl;
			continue;
		}
		testCalcStateIdentificationSets(minFsm, "TC-Rand-(MSP)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(MSP)-" + to_string(i));
	}
	srand(436);
	for (int i = 0; i < 2500; ++i) {
		auto m = Dfsm("M", rand() % 15 + 1, rand() % 4, rand() % 4 + 1, pl, true);
		auto minM = m.minimise();
		testCalcStateIdentificationSets(minM, "TC-Rand-(Dfsm)-" + to_string(i)) 
			? ++result.pass : result.fails.push_back("TC-Rand-(Dfsm)-" + to_string(i));
	}
	srand(63129);
	for (int i = 0; i < 2500; ++i) {
		auto m = Dfsm(*makeStatesPartial(make_shared<Dfsm>(Dfsm("M", rand() % 15 + 1, rand() % 4 + 1, rand() % 4 + 1, pl, true))));
		auto minM = m.minimise();
		testCalcStateIdentificationSets(minM, "TC-Rand-(Dfsm,MSP)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(Dfsm,MSP)-" + to_string(i));
	}

	cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
	Dfsm csm = Dfsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
	testCalcStateIdentificationSets(csm, "TC-CSM-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");

	cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
	Dfsm fsb = Dfsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
	testCalcStateIdentificationSets(fsb, "TC-FSBC-0") ? ++result.pass : result.fails.push_back("TC-FSBC-0");

	cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
	Dfsm gdc = Dfsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
	testCalcStateIdentificationSets(gdc, "TC-GDC-0") ? ++result.pass : result.fails.push_back("TC-GDC-0");

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	vector<DistinguishingTraceTestCase<Fsm>> testSuite;
	parseDistinguishingTraceTSFile("../../../resources/TestSuites/DistinguishingTraces/Fsm_calcStateIdentificationSets.testsuite", testSuite);
	for (auto tc : testSuite) {
		testCalcStateIdentificationSets(*tc.m, "TC-Part-" + tc.id) ? ++result.pass : result.fails.push_back("TC-Part-" + tc.id);
	}
	return result;
}


/**
 * Create a complete Fsm from m by adding self loops in states for undefined inputs producing some nullouput not contained in the
 * regular output alphabet. (nullOutput is expected to be greater than m.getMaxOutput())
 */
shared_ptr<Fsm> transformToComplete(const shared_ptr<const Fsm> m, const size_t nullOutput) {
	vector<shared_ptr<FsmNode> > lst;
	for (int n = 0; n <= m->getMaxState(); n++) {
		lst.push_back(make_shared<FsmNode>(n, m->getName(), m->getPresentationLayer()));
	}

	bool partial = false;

	// Now add transitions that correspond exactly to the transitions in
	// this FSM
	for (int n = 0; n <= m->getMaxState(); n++) {
		auto theNewFsmNodeSrc = lst[n];
		auto theOldFsmNodeSrc = m->getNodes()[n];
		set<int> definedInputs;
		for (auto tr : theOldFsmNodeSrc->getTransitions()) {
			definedInputs.insert(tr->getLabel()->getInput());
			int tgtId = tr->getTarget()->getId();
			auto newLbl = make_shared<FsmLabel>(*(tr->getLabel()));
			shared_ptr<FsmTransition> newTr =
				make_shared<FsmTransition>(theNewFsmNodeSrc, lst[tgtId], newLbl);
			theNewFsmNodeSrc->addTransition(newTr);
		}
		// add self loops with nullOutputs for undefined inputs
		for (int input = 0; input <= m->getMaxInput(); ++input) {
			if (definedInputs.count(input) == 0) {
				partial = true;
				shared_ptr<FsmTransition> newTr =
					make_shared<FsmTransition>(theNewFsmNodeSrc, lst[n], make_shared<FsmLabel>(input, nullOutput, m->getPresentationLayer()));
				theNewFsmNodeSrc->addTransition(newTr);
			}
		}
	}
	if (partial) {
		auto completeM = make_shared<Fsm>(m->getName(), m->getMaxInput(), nullOutput, lst, m->getPresentationLayer());
		completeM->initStateIdx = m->initStateIdx;
		return completeM;
	}
	else {
		auto completeM = make_shared<Fsm>(m->getName(), m->getMaxInput(), m->getMaxOutput(), lst, m->getPresentationLayer());
		completeM->initStateIdx = m->initStateIdx;
		return completeM;
	}
}

/*
 *	Random Test Suite for test of Fsm::calcStateIdentificationSetsFast().
*/
TestResult calcStateIdentificationSetsFast_TS() {
	TestResult result("Fsm::calcStateIdentificationSetsFast");
	cout << "============================= Start Test of Fsm::calcStateIdentificationSetsFast =============================" << endl;

	cout << "------------------------------- Start Random Tests -------------------------------" << endl;	
	srand(1376);
	for (int i = 0; i < 2500; ++i) {
		auto fsm = Fsm::createRandomFsmRepeatable("M1", rand() % 4, rand() % 4 + 1, rand() % 6, make_shared<FsmPresentationLayer>());
		auto minFsm = fsm->minimise();
		if (minFsm.size() > 50) {
			cout << "M is too big. Stop Test Case." << endl;
			continue;
		}
		testCalcStateIdentificationSetsFast(minFsm, "TC-Rand-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}
	srand(711100);
	for (int i = 0; i < 2500; ++i) {
		auto fsm = Fsm::createRandomFsmRepeatable("M1", rand() % 4 + 1, rand() % 4 + 1, rand() % 6, make_shared<FsmPresentationLayer>());
		auto tmp = makeStatesPartial(fsm);
		auto minFsm = tmp->minimise();
		if (minFsm.size() > 50) {
			cout << "M is too big. Stop Test Case." << endl;
			continue;
		}
		testCalcStateIdentificationSetsFast(minFsm, "TC-Rand-(MSP)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(MSP)-" + to_string(i));
	}
	srand(266306);
	for (int i = 0; i < 2500; ++i) {
		auto m = Dfsm("M", rand() % 15 + 1, rand() % 4, rand() % 4 + 1, make_shared<FsmPresentationLayer>(), true);
		auto minM = m.minimise();
		testCalcStateIdentificationSetsFast(minM, "TC-Rand-(Dfsm)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(Dfsm)-" + to_string(i));
	}
	srand(485639);
	for (int i = 0; i < 2500; ++i) {
		auto m = Dfsm("M", rand() % 15 + 1, rand() % 4 + 1, rand() % 4 + 1, make_shared<FsmPresentationLayer>(), true);
		Dfsm minM = Dfsm(*makeStatesPartial(make_shared<Dfsm>(m))).minimise();
		testCalcStateIdentificationSetsFast(minM, "TC-Rand-(Dfsm,MSP)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(Dfsm,MSP)-" + to_string(i));
	}

	cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
	Dfsm csm = Dfsm("../../../resources/TestSuites/examples/csm.fsm", make_shared<FsmPresentationLayer>(), "CSM").minimise();
	testCalcStateIdentificationSetsFast(csm, "TC-CSM-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");

	cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
	Dfsm fsb = Dfsm("../../../resources/TestSuites/examples/fsb.fsm", make_shared<FsmPresentationLayer>(), "FSB").minimise();
	testCalcStateIdentificationSetsFast(fsb, "TC-FSBC-0") ? ++result.pass : result.fails.push_back("TC-FSBC-0");

	cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
	Dfsm gdc = Dfsm("../../../resources/TestSuites/examples/gdc.fsm", make_shared<FsmPresentationLayer>(), "GDC").minimise();
	testCalcStateIdentificationSetsFast(gdc, "TC-GDC-0") ? ++result.pass : result.fails.push_back("TC-GDC-0");


	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	vector<DistinguishingTraceTestCase<Fsm>> testSuite;
	parseDistinguishingTraceTSFile("../../../resources/TestSuites/DistinguishingTraces/Fsm_calcStateIdentificationSetsFast.testsuite", testSuite);
	for (auto tc : testSuite) {
		testCalcStateIdentificationSetsFast(*tc.m, "TC-Part-" + tc.id) ? ++result.pass : result.fails.push_back("TC-Part-" + tc.id);
	}
	return result;
}

// =====================================================================================================================================

// Test functions defined for implemented testing theories.

/**
 * Calculate and return the maxOutput of the Fsms m and mutants.
 */
int getMaxOutput(const Fsm & m, const vector<shared_ptr<const Fsm>>& mutants) {
	int maxOutput = m.getMaxOutput();
	for (shared_ptr<const Fsm > mut : mutants) {
		maxOutput < mut->getMaxOutput() ? maxOutput = mut->getMaxOutput() : maxOutput;
	}
	return maxOutput;
}

class TestSuiteGenerator {
public:
	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates) = 0;
	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates) = 0;
};

class WMethodGenerator : public TestSuiteGenerator {
public:
	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates)
	{
		return make_shared<IOListContainer>(m.wMethod(numAddStates));
	}
	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates)
	{
		return make_shared<IOListContainer>(m.wMethod(numAddStates));
	}
};

class WMethodOnMinimisedFsmGenerator : public TestSuiteGenerator {
public:
	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates)
	{
		return make_shared<IOListContainer>(m.wMethodOnMinimisedFsm(numAddStates));
	}
	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates)
	{
		return make_shared<IOListContainer>(m.wMethodOnMinimisedFsm(numAddStates));
	}
};

class WMethodOnMinimisedDfsmGenerator : public TestSuiteGenerator {
public:
	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates)
	{
		return nullptr;
	}
	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates)
	{
		return make_shared<IOListContainer>(m.wMethodOnMinimisedDfsm(numAddStates));
	}
};

class WpMethodGenerator : public TestSuiteGenerator {
public:
	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates)
	{
		return make_shared<IOListContainer>(m.wpMethod(numAddStates));
	}
	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates)
	{
		return make_shared<IOListContainer>(m.wpMethod(numAddStates));
	}
};

class WpMethodOnMinimisedDfsmGenerator : public TestSuiteGenerator {
public:
	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates)
	{
		return nullptr;
	}
	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates)
	{
		return make_shared<IOListContainer>(m.wpMethodOnMinimisedDfsm(numAddStates));
	}
};

class HsiMethodGenerator : public TestSuiteGenerator {
public:
	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates)
	{
		return make_shared<IOListContainer>(m.hsiMethod(numAddStates));
	}
	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates)
	{
		return make_shared<IOListContainer>(m.hsiMethod(numAddStates));
	}
};

class HMethodOnMinimisedDfsmGenerator : public TestSuiteGenerator {
public:
	virtual shared_ptr<IOListContainer> generateTestSuite(Fsm &m, const unsigned int numAddStates)
	{
		return nullptr;
	}
	virtual shared_ptr<IOListContainer> generateTestSuite(Dfsm &m, const unsigned int numAddStates)
	{
		return make_shared<IOListContainer>(m.hMethodOnMinimisedDfsm(numAddStates));
	}
};

// Calculates and returns numAddStates and adds each mutant that has at most maxAddStates many additional states to filteredMutants
int calcNumAddStatesAndFilterMutants(const vector<shared_ptr<const Fsm>>& mutants, const Fsm& minComplM, vector<shared_ptr<const Fsm>> &filteredMutants) {
	const int maxAddStates = 2;
	int numAddStates = 0;

	for (const auto mutant : mutants) {
		if (mutant->size() > minComplM.size()) {
			int sizeDiff = mutant->size() - minComplM.size();
			if (sizeDiff > maxAddStates) continue;
			filteredMutants.push_back(mutant);
			if (sizeDiff > numAddStates) {
				numAddStates = sizeDiff;
			}
		}
		else {
			filteredMutants.push_back(mutant);
		}
	}
	return numAddStates;
}

// Calculates and compares outputs generated by completeM and mutant for the given test suite. Return true iff there is some test case producing different
// output traces.
bool compareOutputs(const shared_ptr<IOListContainer> ts, const shared_ptr<Fsm> completeM, const shared_ptr<const Fsm> mutant) {
	bool diff = false;
	for (const auto &tc : *ts->getIOLists()) {
		if (calcCompleteOutputTraces(completeM->getInitialState(), tc)
			!= calcCompleteOutputTraces(mutant->getInitialState(), tc)) {
			diff = true;
			break;
		}
	}
	return diff;
}

/**
 * Test function for test theories.
 * Each element of mutants is expected to be complete and normalized and to have the same input alphabet as m.
 * nullOutput is expected to be the output that was used as nulloutputs in the completion of the mutants.
 */
bool testTestTheory(Fsm & m, const vector<shared_ptr<const Fsm>>& mutants, const size_t nullOutput,
	                const shared_ptr<TestSuiteGenerator> tsGen, const string &tcID) {
	bool pass = true;
	const Fsm copyOfM = Fsm(m);		
	auto completeM = transformToComplete(make_shared<Fsm>(m), nullOutput);
	auto minComplM = completeM->minimise();
	if (minComplM.size() > 30) {
		cout << "FSM too big. Stop Test Case." << endl;
		return pass;
	}

	// calculate numAddStates 
	vector<shared_ptr<const Fsm>> filteredMutants;
	int numAddStates = calcNumAddStatesAndFilterMutants(mutants, minComplM, filteredMutants);

	// Calculate complete test suite
	const auto ts = tsGen->generateTestSuite(m, numAddStates);

	// first check invariant of m
	bool invariantViolation = not m.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Invariant still holds for M after calculation.");
	// stop test execution at this point if invariant of m does not hold anymore
	if (invariantViolation) return pass;

	size_t maxLength = 0;
	for (auto tc : *ts->getIOLists()) { if (tc.size() > maxLength) maxLength = tc.size(); }

	// stop test case if test cases are too long or test suite contains too many test cases
	if (maxLength > 12 or ts->size() > 2000) {
		cout << "Test Suite too big or test cases too long. Stop test case." << endl;
		return pass;
	}

	// Check completeness of test suite with help of the mutants
	for (const auto mutant : filteredMutants) {
		bool diff = compareOutputs(ts, completeM, mutant);
		bool eq = ioEquivalenceCheck(completeM->getInitialState(), mutant->getInitialState());
		fsmlib_assert(tcID, eq != diff, pass, "M and mutant are i/o-equivalent iff mutant passed test suite.");
		if (eq == diff) {
			cout << *completeM << endl;
			cout << *mutant << endl;
			cout << "TS: " << *ts << endl;
			cout << "nullOutput: " << nullOutput << endl;
		}
	}

	// check if language of m has changed
	fsmlib_assert(tcID, ioEquivalenceCheck(m.getInitialState(), copyOfM.getInitialState()), pass, "Language of M has not changed");
	return pass;
}

/**
 * Test function for Fsm::wMethod(const unsigned int numAddStates)
 * Each element of mutants is expected to be complete and normalized and to have the same input alphabet as m.
 * nullOutput is expected to be the output that was used as nulloutputs in the completion of the mutants.
 */
bool testTestTheory(Dfsm & m, const vector<shared_ptr<const Fsm>>& mutants, const size_t nullOutput, 
	                const shared_ptr<TestSuiteGenerator> tsGen, const string &tcID) {
	bool pass = true;
	const Dfsm copyOfM = Dfsm(m);
	auto completeM = transformToComplete(make_shared<Dfsm>(m), nullOutput);
	auto minComplM = completeM->minimise();
	if (minComplM.size() > 30) {
		cout << "FSM too big. Stop Test Case." << endl;
		return pass;
	}
	// calculate numAddStates 
	vector<shared_ptr<const Fsm>> filteredMutants;
	int numAddStates = calcNumAddStatesAndFilterMutants(mutants, minComplM, filteredMutants);

	// Calculate complete test suite
	const auto ts = tsGen->generateTestSuite(m, numAddStates);

	// first check invariant of m
	bool invariantViolation = not m.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Invariant still holds for M after calculation.");
	// stop test execution at this point if invariant of m does not hold anymore
	if (invariantViolation) return pass;

	size_t maxLength = 0;
	for (auto tc : *ts->getIOLists()) { if (tc.size() > maxLength) maxLength = tc.size(); }

	for (const auto mutant : filteredMutants) {
		bool diff = compareOutputs(ts, completeM, mutant);
		bool eq = ioEquivalenceCheck(completeM->getInitialState(), mutant->getInitialState());
		fsmlib_assert(tcID, eq != diff, pass, "M and mutant are i/o-equivalent iff mutant passed test suite.");
		if (eq == diff) {
			cout << *completeM << endl;
			cout << *mutant << endl;
			cout << "TS: " << *ts << endl;
			cout << "nullOutput: " << nullOutput << endl;
		}
	}

	// check if language of m has changed
	fsmlib_assert(tcID, ioEquivalenceCheck(m.getInitialState(), copyOfM.getInitialState()), pass, "Language of M has not changed");
	return pass;
}

struct TestTheoryTestCase {
	string id;
	string rmPath;
	vector<string> mutantPaths;
};

shared_ptr<vector<TestTheoryTestCase>> parseTestTheoryTSFile(const string &testSuitePath) {
	string fname = testSuitePath;
	vector<TestTheoryTestCase> testSuite;
	ifstream inputFile(fname);
	if (inputFile.is_open())
	{
		string line;
		while (getline(inputFile, line))
		{
			vector<string> lineContent;
			stringstream ss(line);
			for (string elem; getline(ss, elem, ';'); lineContent.push_back(elem));
			TestTheoryTestCase tc;
			tc.id = lineContent.at(0);
			tc.rmPath = lineContent.at(1);
			for (size_t i = 2; i < lineContent.size(); ++i) {
				tc.mutantPaths.push_back(lineContent.at(i));
			}
			testSuite.push_back(tc);
		}
		inputFile.close();
		return make_shared<vector<TestTheoryTestCase>>(testSuite);
	}
	else
	{
		cout << "Unable to open input file" << endl;
		exit(EXIT_FAILURE);
	}
}

// Creates minimised and completed mutants from m.
shared_ptr<vector<shared_ptr<const Fsm>>> createMutants(const size_t nullOutput, const shared_ptr<Fsm> m) {
	vector<shared_ptr<const Fsm>> mutants;
	for (int j = 0; j < 20; ++j) {
		size_t numOutFaults = rand() % 3;
		size_t numTrFaults = rand() % 3;
		if (numOutFaults == 0 and numTrFaults == 0) ++numTrFaults;  // ignore the case where both values equal 0
		auto minMut = m->createMutantRepeatable("Mutant_" + to_string(j), numOutFaults, numTrFaults)->minimise();
		mutants.push_back(transformToComplete(make_shared<Fsm>(minMut), nullOutput));
	}
	return make_shared< vector<shared_ptr<const Fsm>>>(mutants);
}

template<typename T>
void executeTestTheoryTC(TestTheoryTestCase & tc, shared_ptr<TestSuiteGenerator> tsGenerator, TestResult &result) {
	shared_ptr<T> ref = make_shared<T>(tc.rmPath, make_shared<FsmPresentationLayer>(), "M");
	vector<shared_ptr<const Fsm>> partialMutants;
	for (size_t i = 0; i < tc.mutantPaths.size(); ++i) {
		shared_ptr<const Fsm> mutant = make_shared<Fsm>(tc.mutantPaths.at(i), make_shared<FsmPresentationLayer>(), "Mut_" + to_string(i));
		partialMutants.push_back(mutant);
	}
	size_t nullOutput = getMaxOutput(*ref, partialMutants) + 1;

	vector<shared_ptr<const Fsm>> completeMutants;
	for (auto m : partialMutants) completeMutants.push_back(transformToComplete(m, nullOutput));
	testTestTheory(*ref, completeMutants, nullOutput, tsGenerator, "TC-Part-" + tc.id) ? ++result.pass : result.fails.push_back("TC-Part-" + tc.id);
}


// Test Fsm::wMethod(...)
TestResult wMethod_Fsm_TS() {
	TestResult result("Fsm::wMethod");
	cout << "============================= Start Test of Fsm::wMethod =============================" << endl;
	
	shared_ptr<WMethodGenerator> tsGenerator = make_shared<WMethodGenerator>();
	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	srand(12792);
	for (int i = 0; i < 2000; ++i) {
		auto m = makeStatesUnreachable(*Fsm::createRandomFsmRepeatable("M", (rand() % 4) + 1, (rand() % 6) + 1, (rand() % 6) + 1, pl));
		const size_t nullOutput = m->getMaxOutput() + 1;

		testTestTheory(*m, *createMutants(nullOutput, m), nullOutput, tsGenerator, "TC-Rand-(MSU)-"+ to_string(i)) 
			? ++result.pass : result.fails.push_back("TC-Rand-(MSU)-" + to_string(i));
	}

	srand(94563);
	for (int i = 0; i < 2000; ++i) {
		auto m = Fsm::createRandomFsmRepeatable("M", (rand() % 4) + 1, (rand() % 6) + 1, (rand() % 6) + 1, pl);
		const size_t nullOutput = m->getMaxOutput() + 1;

		testTestTheory(*m, *createMutants(nullOutput, m), nullOutput, tsGenerator, "TC-Rand-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}
	srand(3641);
	for (int i = 0; i < 2000; ++i) {
		auto m = makeStatesEquivalent(*Fsm::createRandomFsmRepeatable("M", (rand() % 4) + 1, (rand() % 6) + 1, (rand() % 6) + 1, pl));
		const size_t nullOutput = m->getMaxOutput() + 1;

		testTestTheory(*m, *createMutants(nullOutput, m), nullOutput, tsGenerator, "TC-Rand-(MSE)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(MSE)-" + to_string(i));
	}
	srand(94775);
	for (int i = 0; i < 2000; ++i) {
		auto m = makeStatesPartial(Fsm::createRandomFsmRepeatable("M", (rand() % 4) + 1, (rand() % 6) + 1, (rand() % 6) + 1, pl));
		const size_t nullOutput = m->getMaxOutput() + 1;

		testTestTheory(*m, *createMutants(nullOutput, m), nullOutput, tsGenerator, "TC-Rand-(MSP)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(MSP)-" + to_string(i));
	}
	srand(2022);
	for (int i = 0; i < 2000; ++i) {
		auto m = makeStatesEquivalent(*makeStatesPartial(Fsm::createRandomFsmRepeatable("M", (rand() % 4) + 1, (rand() % 6) + 1, (rand() % 6) + 1, pl)));
		const size_t nullOutput = m->getMaxOutput() + 1;

		testTestTheory(*m, *createMutants(nullOutput, m), nullOutput, tsGenerator, "TC-Rand-(MSE,MSP)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(MSE,MSP)-" + to_string(i));
	}
	srand(52151);
	for (int i = 0; i < 2000; ++i) {
		auto m = makeStatesPartial(makeStatesUnreachable(*Fsm::createRandomFsmRepeatable("M", (rand() % 4) + 1, (rand() % 6) + 1, (rand() % 6) + 1, pl)));
		const size_t nullOutput = m->getMaxOutput() + 1;

		testTestTheory(*m, *createMutants(nullOutput, m), nullOutput, tsGenerator, "TC-Rand-(MSU)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(MSU,MSP)-" + to_string(i));
	}

	srand(990137);
	{
		cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
		shared_ptr<Fsm> csm = make_shared<Fsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// CSM is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(csm->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(*csm, mutants, csm->getMaxOutput() + 1, tsGenerator, "TC-CSM-0") 
			? ++result.pass : result.fails.push_back("TC-CSM-0");
	}
	{
		cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
		shared_ptr<Fsm> fsbc = make_shared<Fsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// FSB is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(fsbc->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(*fsbc, mutants, fsbc->getMaxOutput() + 1, tsGenerator, "TC-FSBC-0")
			? ++result.pass : result.fails.push_back("TC-FSBC-0");
	}
	{
		cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
		shared_ptr<Fsm> gdc = make_shared<Fsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// GDC is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(gdc->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(*gdc, mutants, gdc->getMaxOutput() + 1, tsGenerator, "TC-GDC-0")
			? ++result.pass : result.fails.push_back("TC-GDC-0");
	}

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Fsm_wMethod.testsuite");
	for (auto tc : *testSuite) {
		executeTestTheoryTC<Fsm>(tc, tsGenerator, result);
	}
	return result;
}

// Test Dfsm::wMethod(...)
TestResult wMethod_Dfsm_TS() {
	TestResult result("Dfsm::wMethod");
	cout << "============================= Start Test of Dfsm::wMethod =============================" << endl;	
	shared_ptr<WMethodGenerator> tsGenerator = make_shared<WMethodGenerator>();
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(412725);
	for (int i = 0; i < 100; ++i) {
		Dfsm m("M", rand() % 15 + 1, rand() % 6, (rand() % 6) + 1, pl, true);
		const size_t nullOutput = m.getMaxOutput() + 1;
		testTestTheory(m, *createMutants(nullOutput, make_shared<Dfsm>(m)), nullOutput, tsGenerator, "TC-Rand-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}

	srand(130343);
	{
		cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
		shared_ptr<Dfsm> csm = make_shared<Dfsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// CSM is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(csm->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(*csm, mutants, csm->getMaxOutput() + 1, tsGenerator, "TC-CSM-0")
			? ++result.pass : result.fails.push_back("TC-CSM-0");
	}
	{
		cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
		shared_ptr<Dfsm> fsb = make_shared<Dfsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// FSB is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(fsb->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(*fsb, mutants, fsb->getMaxOutput() + 1, tsGenerator, "TC-FSBC-0")
			? ++result.pass : result.fails.push_back("TC-FSBC-0");
	}
	{
		cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
		shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// GDC is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(gdc->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(*gdc, mutants, gdc->getMaxOutput() + 1, tsGenerator, "TC-GDC-0")
			? ++result.pass : result.fails.push_back("TC-GDC-0");
	}

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Dfsm_wMethod.testsuite");
	for (auto tc : *testSuite) {
		executeTestTheoryTC<Dfsm>(tc, tsGenerator, result);
	}
	return result;
}

// Test Fsm::wMethodOnMinimisedFsm(...)
TestResult wMethodOnMinimisedFsm_TS() {
	TestResult result("Fsm::wMethodOnMinimisedFsm");
	cout << "============================= Start Test of Fsm::wMethodOnMinimisedFsm =============================" << endl;	
	shared_ptr<WMethodOnMinimisedFsmGenerator> tsGenerator = make_shared<WMethodOnMinimisedFsmGenerator>();
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(50301);
	for (int i = 0; i < 2500; ++i) {
		auto m = Fsm::createRandomFsmRepeatable("M", rand() % 6, (rand() % 6) + 1, rand() % 6 + 1, pl)->minimise();
		// filter fsm that have too many states
		if (m.size() > 20) {
			cout << "FSM too big. Stop Test Case." << endl;
			continue;
		}		
		const size_t nullOutput = m.getMaxOutput() + 1;
		testTestTheory(m, *createMutants(nullOutput, make_shared<Fsm>(m)), nullOutput, tsGenerator, "TC-Rand-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}
	srand(27161);
	for (int i = 0; i < 2500; ++i) {
		auto m = makeStatesPartial(Fsm::createRandomFsmRepeatable("M", rand() % 6 + 1, (rand() % 6) + 1, rand() % 6 + 1, pl))->minimise();
		// filter fsm that have too many states
		if (m.size() > 20) {
			cout << "FSM too big. Stop Test Case." << endl;
			continue;
		}
		const size_t nullOutput = m.getMaxOutput() + 1;
		testTestTheory(m, *createMutants(nullOutput, make_shared<Fsm>(m)), nullOutput, tsGenerator, "TC-Rand-(MSP)-" + to_string(i)) 
			? ++result.pass : result.fails.push_back("TC-Rand-(MSP)-" + to_string(i));
	}
	srand(775565);
	for (int i = 0; i < 2500; ++i) {
		Dfsm m("M", rand() % 15 + 1, rand() % 6, (rand() % 6) + 1, pl, true);
		m = m.minimise();

		const size_t nullOutput = m.getMaxOutput() + 1;
		testTestTheory(m, *createMutants(nullOutput, make_shared<Dfsm>(m)), nullOutput, tsGenerator, "TC-Rand-(Dfsm)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(Dfsm)-" + to_string(i));
	}
	srand(67131);
	for (int i = 0; i < 2500; ++i) {
		Dfsm m(*makeStatesPartial(make_shared<Dfsm>(Dfsm("M", rand() % 15 + 1, rand() % 6 + 1, (rand() % 6) + 1, pl, true))));
		m = m.minimise();

		const size_t nullOutput = m.getMaxOutput() + 1;
		testTestTheory(m, *createMutants(nullOutput, make_shared<Dfsm>(m)), nullOutput, tsGenerator, "TC-Rand-(Dfsm,MSP)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(Dfsm,MSP)-" + to_string(i));
	}

	srand(326744);
	{
		cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
		Dfsm csm = Dfsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// CSM is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(csm.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(csm, mutants, csm.getMaxOutput() + 1, tsGenerator, "TC-CSM-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");
	}
	{
		cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
		Dfsm fsb = Dfsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// FSB is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(fsb.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(fsb, mutants, fsb.getMaxOutput() + 1, tsGenerator, "TC-FSBC-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");
	}
	{
		cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
		Dfsm gdc = Dfsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// GDC is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(gdc.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(gdc, mutants, gdc.getMaxOutput() + 1, tsGenerator, "TC-GDC-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");
	}

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Fsm_wMethodOnMinimisedFsm.testsuite");
	for (auto tc : *testSuite) {
		executeTestTheoryTC<Fsm>(tc, tsGenerator, result);
	}
	return result;
}

// test Dfsm::wMethodOnMinimisedDfsm(...)
TestResult wMethodOnMinimisedDfsm_TS() {
	TestResult result("Dfsm::wMethodOnMinimisedDfsm");
	cout << "============================= Start Test of Dfsm::wMethodOnMinimisedDfsm =============================" << endl;	
	shared_ptr<WMethodOnMinimisedDfsmGenerator> tsGenerator = make_shared<WMethodOnMinimisedDfsmGenerator>();
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(625);
	for (int i = 0; i < 10000; ++i) {
		Dfsm m("M", rand() % 15 + 1, rand() % 6, (rand() % 6) + 1, pl, true);
		m = m.minimise();
		const size_t nullOutput = m.getMaxOutput() + 1;
		testTestTheory(m, *createMutants(nullOutput, make_shared<Dfsm>(m)), nullOutput, tsGenerator, "TC-Rand-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}
	srand(69080);
	{
		cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
		Dfsm csm = Dfsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// CSM is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(csm.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(csm, mutants, csm.getMaxOutput() + 1, tsGenerator, "TC-CSM-0")
			? ++result.pass : result.fails.push_back("TC-CSM-0");
	}
	{
		cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
		Dfsm fsb = Dfsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// FSB is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(fsb.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(fsb, mutants, fsb.getMaxOutput() + 1, tsGenerator, "TC-FSBC-0")
			? ++result.pass : result.fails.push_back("TC-FSBC-0");
	}
	{
		cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
		Dfsm gdc = Dfsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// GDC is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(gdc.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(gdc, mutants, gdc.getMaxOutput() + 1, tsGenerator, "TC-GDC-0")
			? ++result.pass : result.fails.push_back("TC-GDC-0");
	}

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Dfsm_wMethodOnMinimisedDfsm.testsuite");
	for (auto tc : *testSuite) {
		executeTestTheoryTC<Dfsm>(tc, tsGenerator, result);
	}
	return result;
}

// test Fsm::wpMethod(...)
TestResult wpMethod_Fsm_TS() {
	TestResult result("Fsm::wpMethod");
	cout << "============================= Start Test of Fsm::wpMethod =============================" << endl;	
	shared_ptr<WpMethodGenerator> tsGenerator = make_shared<WpMethodGenerator>();
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(5217);
	for (int i = 0; i < 5000; ++i) {
		auto m = Fsm::createRandomFsmRepeatable("M", rand() % 4, (rand() % 6) + 1, rand() % 6 + 1, pl)->minimise();
		const size_t nullOutput = m.getMaxOutput() + 1;

		// filter fsm that have too many states
		if (m.size() > 20) {
			cout << "FSM too big. Stop Test Case." << endl;
			continue;
		}
		testTestTheory(m, *createMutants(nullOutput, make_shared<Fsm>(m)), nullOutput, tsGenerator, "TC-Rand-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}
	srand(41575);
	for (int i = 0; i < 5000; ++i) {
		auto m = makeStatesPartial(Fsm::createRandomFsmRepeatable("M", rand() % 4 + 1, (rand() % 6) + 1, rand() % 6 + 1, pl))->minimise();
		const size_t nullOutput = m.getMaxOutput() + 1;

		// filter fsm that have too many states
		if (m.size() > 20) {
			cout << "FSM too big. Stop Test Case." << endl;
			continue;
		}
		testTestTheory(m, *createMutants(nullOutput, make_shared<Fsm>(m)), nullOutput, tsGenerator, "TC-Rand-(MSP)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(MSP)-" + to_string(i));
	}

	srand(156041);
	{
		cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
		Fsm csm = Fsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// CSM is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(csm.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(csm, mutants, csm.getMaxOutput() + 1, tsGenerator, "TC-CSM-0")
			? ++result.pass : result.fails.push_back("TC-CSM-0");
	}
	{
		cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
		Fsm fsb = Fsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// FSB is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(fsb.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(fsb, mutants, fsb.getMaxOutput() + 1, tsGenerator, "TC-FSBC-0")
			? ++result.pass : result.fails.push_back("TC-FSBC-0");
	}
	{
		cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
		Fsm gdc = Fsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// GDC is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(gdc.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(gdc, mutants, gdc.getMaxOutput() + 1, tsGenerator, "TC-GDC-0")
			? ++result.pass : result.fails.push_back("TC-GDC-0");
	}

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Fsm_wpMethod.testsuite");
	for (auto tc : *testSuite) {
		executeTestTheoryTC<Fsm>(tc, tsGenerator, result);
	}
	return result;
}

// Test Dfsm::wpMethod(...)
TestResult wpMethod_Dfsm_TS() {
	TestResult result("Dfsm::wpMethod");
	cout << "============================= Start Test of Dfsm::wpMethod =============================" << endl;	
	shared_ptr<WpMethodGenerator> tsGenerator = make_shared<WpMethodGenerator>();
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(712871);
	for (int i = 0; i < 10000; ++i) {
		Dfsm m("M", rand() % 15 + 1, rand() % 6, (rand() % 6) + 1, pl, true);

		const size_t nullOutput = m.getMaxOutput() + 1;

		testTestTheory(m, *createMutants(nullOutput, make_shared<Dfsm>(m)), nullOutput, tsGenerator, "TC-Rand-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}

	srand(37667);
	{
		cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
		shared_ptr<Dfsm> csm = make_shared<Dfsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// CSM is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(csm->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(*csm, mutants, csm->getMaxOutput() + 1, tsGenerator, "TC-CSM-0")
			? ++result.pass : result.fails.push_back("TC-CSM-0");
	}
	{
		cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
		shared_ptr<Dfsm> fsb = make_shared<Dfsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// FSB is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(fsb->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(*fsb, mutants, fsb->getMaxOutput() + 1, tsGenerator, "TC-FSBC-0")
			? ++result.pass : result.fails.push_back("TC-FSBC-0");
	}
	{
		cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
		shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// GDC is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(gdc->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(*gdc, mutants, gdc->getMaxOutput() + 1, tsGenerator, "TC-GDC-0")
			? ++result.pass : result.fails.push_back("TC-GDC-0");
	}

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Dfsm_wpMethod.testsuite");
	for (auto tc : *testSuite) {
		executeTestTheoryTC<Dfsm>(tc, tsGenerator, result);
	}
	return result;
}

// Test Dfsm::wpMethodOnMinimisedDfsm(...)
TestResult wpMethodOnMinimisedDfsm_TS() {
	TestResult result("Dfsm::wpMethodOnMinimisedDfsm");
	cout << "============================= Start Test of Dfsm::wpMethodOnMinimisedDfsm =============================" << endl;	
	shared_ptr<WpMethodOnMinimisedDfsmGenerator> tsGenerator = make_shared<WpMethodOnMinimisedDfsmGenerator>();
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(625166);
	for (int i = 0; i < 10000; ++i) {
		Dfsm m("M", rand() % 15 + 1, rand() % 6, (rand() % 6) + 1, pl, true);
		m = m.minimise();

		const size_t nullOutput = m.getMaxOutput() + 1;

		testTestTheory(m, *createMutants(nullOutput, make_shared<Dfsm>(m)), nullOutput, tsGenerator, "TC-Rand-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}

	srand(81460);
	{
		cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
		Dfsm csm = Dfsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// CSM is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(csm.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(csm, mutants, csm.getMaxOutput() + 1, tsGenerator, "TC-CSM-0")
			? ++result.pass : result.fails.push_back("TC-CSM-0");
	}
	{
		cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
		Dfsm fsb = Dfsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// FSB is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(fsb.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(fsb, mutants, fsb.getMaxOutput() + 1, tsGenerator, "TC-FSBC-0")
			? ++result.pass : result.fails.push_back("TC-FSBC-0");
	}
	{
		cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
		Dfsm gdc = Dfsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// GDC is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(gdc.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(gdc, mutants, gdc.getMaxOutput() + 1, tsGenerator, "TC-GDC-0")
			? ++result.pass : result.fails.push_back("TC-GDC-0");
	}

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Dfsm_wpMethodOnMinimisedDfsm.testsuite");
	for (auto tc : *testSuite) {
		executeTestTheoryTC<Dfsm>(tc, tsGenerator, result);
	}
	return result;
}

// test Fsm::hsiMethod(...)
TestResult hsiMethod_Fsm_TS() {
	TestResult result("Fsm::hsiMethod");
	cout << "============================= Start Test of Fsm::hsiMethod =============================" << endl;	
	shared_ptr<HsiMethodGenerator> tsGenerator = make_shared<HsiMethodGenerator>();
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(42676);
	for (int i = 0; i < 5000; ++i) {
		auto m = Fsm::createRandomFsmRepeatable("M", rand() % 4, (rand() % 6) + 1, rand() % 6 + 1, pl)->minimise();
		const size_t nullOutput = m.getMaxOutput() + 1;

		// filter fsm that have too many states
		if (m.size() > 20) {
			cout << "FSM too big. Stop Test Case." << endl;
			continue;
		}
		testTestTheory(m, *createMutants(nullOutput, make_shared<Fsm>(m)), nullOutput, tsGenerator, "TC-Rand-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}
	srand(12);
	for (int i = 0; i < 5000; ++i) {
		auto m = makeStatesPartial(Fsm::createRandomFsmRepeatable("M", rand() % 4 + 1, (rand() % 6) + 1, rand() % 6 + 1, pl))->minimise();
		const size_t nullOutput = m.getMaxOutput() + 1;

		// filter fsm that have too many states
		if (m.size() > 20) {
			cout << "FSM too big. Stop Test Case." << endl;
			continue;
		}
		testTestTheory(m, *createMutants(nullOutput, make_shared<Fsm>(m)), nullOutput, tsGenerator, "TC-Rand-(MSP)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(MSP)-" + to_string(i));
	}

	srand(15937);
	{
		cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
		Fsm csm = Fsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// CSM is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(csm.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(csm, mutants, csm.getMaxOutput() + 1, tsGenerator, "TC-CSM-0")
			? ++result.pass : result.fails.push_back("TC-CSM-0");
	}
	{
		cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
		Fsm fsb = Fsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// FSB is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(fsb.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(fsb, mutants, fsb.getMaxOutput() + 1, tsGenerator, "TC-FSBC-0")
			? ++result.pass : result.fails.push_back("TC-FSBC-0");
	}
	{
		cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
		Fsm gdc = Fsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// GDC is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(gdc.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(gdc, mutants, gdc.getMaxOutput() + 1, tsGenerator, "TC-GDC-0")
			? ++result.pass : result.fails.push_back("TC-GDC-0");
	}

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Fsm_hsiMethod.testsuite");
	for (auto tc : *testSuite) {
		executeTestTheoryTC<Fsm>(tc, tsGenerator, result);
	}
	return result;
}

// Test Dfsm::hsiMethod(...)
TestResult hsiMethod_Dfsm_TS() {
	TestResult result("Dfsm::hsiMethod");
	cout << "============================= Start Test of Dfsm::hsiMethod =============================" << endl;	
	shared_ptr<HsiMethodGenerator> tsGenerator = make_shared<HsiMethodGenerator>();
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(907);
	for (int i = 0; i < 5000; ++i) {
		Dfsm m("M", rand() % 15 + 1, rand() % 6, (rand() % 6) + 1, pl, true);

		const size_t nullOutput = m.getMaxOutput() + 1;

		testTestTheory(m, *createMutants(nullOutput, make_shared<Dfsm>(m)), nullOutput, tsGenerator, "TC-Rand-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}
	srand(6144);
	for (int i = 0; i < 5000; ++i) {
		Dfsm m(*makeStatesPartial(make_shared<Dfsm>(Dfsm("M", rand() % 15 + 1, rand() % 6 + 1, (rand() % 6) + 1, pl, true))));

		const size_t nullOutput = m.getMaxOutput() + 1;

		testTestTheory(m, *createMutants(nullOutput, make_shared<Dfsm>(m)), nullOutput, tsGenerator, "TC-Rand-(MSP)-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-(MSP)-" + to_string(i));
	}

	srand(349);
	{
		cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
		shared_ptr<Dfsm> csm = make_shared<Dfsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// CSM is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(csm->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(*csm, mutants, csm->getMaxOutput() + 1, tsGenerator, "TC-CSM-0")
			? ++result.pass : result.fails.push_back("TC-CSM-0");
	}
	{
		cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
		shared_ptr<Dfsm> fsb = make_shared<Dfsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// FSB is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(fsb->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(*fsb, mutants, fsb->getMaxOutput() + 1, tsGenerator, "TC-FSBC-0")
			? ++result.pass : result.fails.push_back("TC-FSBC-0");
	}
	{
		cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
		shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// GDC is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(gdc->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(*gdc, mutants, gdc->getMaxOutput() + 1, tsGenerator, "TC-GDC-0")
			? ++result.pass : result.fails.push_back("TC-GDC-0");
	}

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Dfsm_hsiMethod.testsuite");
	for (auto tc : *testSuite) {
		executeTestTheoryTC<Dfsm>(tc, tsGenerator, result);
	}
	return result;
}

// Test Dfsm::hMethodOnMinimisedDfsm(...)
TestResult hMethodOnMinimisedDfsm_TS() {
	TestResult result("Dfsm::hMethodOnMinimisedDfsm");
	cout << "============================= Start Test of Dfsm::hMethodOnMinimisedDfsm =============================" << endl;	
	shared_ptr<HMethodOnMinimisedDfsmGenerator> tsGenerator = make_shared<HMethodOnMinimisedDfsmGenerator>();
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(3234);
	for (int i = 0; i < 10000; ++i) {
		Dfsm m("M", rand() % 13 + 2, rand() % 6 + 1, (rand() % 6) + 1, pl, true);
		m = m.minimise();
		if (m.size() < 2) {
			cout << "M has less than 2 states. Stop test case." << endl;
			continue;
		}
		const size_t nullOutput = m.getMaxOutput() + 1;
		testTestTheory(m, *createMutants(nullOutput, make_shared<Dfsm>(m)), nullOutput, tsGenerator, "TC-Rand-" + to_string(i))
			? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}

	srand(807);
	{
		cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
		Dfsm csm = Dfsm("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// CSM is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(csm.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(csm, mutants, csm.getMaxOutput() + 1, tsGenerator, "TC-CSM-0")
			? ++result.pass : result.fails.push_back("TC-CSM-0");
	}
	{
		cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
		Dfsm fsb = Dfsm("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// FSB is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(fsb.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(fsb, mutants, fsb.getMaxOutput() + 1, tsGenerator, "TC-FSBC-0")
			? ++result.pass : result.fails.push_back("TC-FSBC-0");
	}
	{
		cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
		Dfsm gdc = Dfsm("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC").minimise();
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			// GDC is complete => created mutants are complete => minimised mutants are complete
			mutants.push_back(make_shared<Fsm>(gdc.createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, rand() % 6 + 1)->minimise()));
		}
		testTestTheory(gdc, mutants, gdc.getMaxOutput() + 1, tsGenerator, "TC-GDC-0")
			? ++result.pass : result.fails.push_back("TC-GDC-0");
	}

	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Dfsm_hMethodOnMinimisedDfsm.testsuite");
	for (auto tc : *testSuite) {
		executeTestTheoryTC<Dfsm>(tc, tsGenerator, result);	
	}
	return result;
}

/**
 * Test function for Dfsm::tMethod().
 * m is expected to be complete. Each element of mutants is expected to differ from m only by zero or more output faults.
 */
bool testTMethod(Dfsm & m, const vector<shared_ptr<const Fsm>>& mutants, const string &tcID) {
	bool pass = true;
	const Dfsm copyOfM = Dfsm(m);

	const auto ts = m.tMethod();

	//int nullOutput = getMaxOutput(m, mutants) + 1;

	// first check invariant of m
	bool invariantViolation = not m.checkInvariant();
	fsmlib_assert(tcID, not invariantViolation, pass, "Dfsm class invariant still holds for M after calculation.");
	// stop test execution at this point if invariant of m does not hold anymore
	if (invariantViolation) return pass;

	for (const auto mutant : mutants) {
		bool diff = false;
		for (const auto &tc : *ts.getIOLists()) {
			if (calcCompleteOutputTraces(m.getInitialState(), tc)
				!= calcCompleteOutputTraces(mutant->getInitialState(), tc)) {
				diff = true;
				break;
			}
		}
		fsmlib_assert(tcID, ioEquivalenceCheck(m.getInitialState(), mutant->getInitialState()) != diff, pass, "M and mutant are i/o-equivalent iff mutant passed test suite.");
	}

	// check if structure of m has changed
	fsmlib_assert(tcID, checkForEqualStructure(m, copyOfM), pass, "M was not changed by algorithm");
	return pass;
}

/*
 *	Test Suite of Dfsm::tMethod().
 */
TestResult tMethod_TS() {
	TestResult result("Dfsm::tMethod");
	cout << "============================= Start Test of Dfsm::tMethod =============================" << endl;	
	shared_ptr<FsmPresentationLayer> pl = make_shared<FsmPresentationLayer>();
	cout << "------------------------------- Start Random Tests -------------------------------" << endl;
	srand(876298);
	for (int i = 0; i < 10000; ++i) {
		auto m = Dfsm("M", (rand() % 15 + 1), rand() % 6, (rand() % 6) + 1, pl, true);
		vector<shared_ptr<const Fsm>> mutants;
		for (int j = 0; j < 20; ++j) {
			size_t numOutFaults = (rand() % 4) + 1;
			mutants.push_back(m.createMutantRepeatable("Mutant_" + to_string(j), numOutFaults, 0));
		}
		testTMethod(m, mutants, "TC-Rand-" + to_string(i)) ? ++result.pass : result.fails.push_back("TC-Rand-" + to_string(i));
	}

	srand(36);
	{
		cout << "------------------------------- Start CSM Tests -------------------------------" << endl;
		shared_ptr<Dfsm> csm = make_shared<Dfsm>("../../../resources/TestSuites/examples/csm.fsm", pl, "CSM");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			mutants.push_back(csm->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, 0));
		}
		testTMethod(*csm, mutants, "TC-CSM-0") ? ++result.pass : result.fails.push_back("TC-CSM-0");
	}
	{
		cout << "------------------------------- Start FSBC Tests -------------------------------" << endl;
		shared_ptr<Dfsm> fsb = make_shared<Dfsm>("../../../resources/TestSuites/examples/fsb.fsm", pl, "FSB");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			mutants.push_back(fsb->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, 0));
		}
		testTMethod(*fsb, mutants, "TC-FSBC-0") ? ++result.pass : result.fails.push_back("TC-FSBC-0");
	}
	{
		cout << "------------------------------- Start GDC Tests -------------------------------" << endl;
		shared_ptr<Dfsm> gdc = make_shared<Dfsm>("../../../resources/TestSuites/examples/gdc.fsm", pl, "GDC");
		vector<shared_ptr<const Fsm>> mutants;
		for (int i = 0; i < 100; ++i) {
			mutants.push_back(gdc->createMutantRepeatable("Mutant_" + to_string(i), rand() % 6 + 1, 0));
		}
		testTMethod(*gdc, mutants, "TC-GDC-0") ? ++result.pass : result.fails.push_back("TC-GDC-0");
	}


	cout << "------------------------------- Start Partition Tests -------------------------------" << endl;
	auto testSuite = parseTestTheoryTSFile("../../../resources/TestSuites/TestTheories/Dfsm_tMethod.testsuite");
	for (auto tc : *testSuite) {
		shared_ptr<Dfsm> ref = make_shared<Dfsm>(tc.rmPath, pl, "M");
		vector<shared_ptr<const Fsm>> completeMutants;
		for (size_t i = 0; i < tc.mutantPaths.size(); ++i) {
			shared_ptr<const Fsm> mutant = make_shared<Fsm>(tc.mutantPaths.at(i), pl, "Mut_" + to_string(i));
			completeMutants.push_back(mutant);
		}
		testTMethod(*ref, completeMutants, "TC-Part-" + tc.id) ? ++result.pass : result.fails.push_back("TC-Part-" + tc.id);
	}
	return result;
}

// ====================================================================================================