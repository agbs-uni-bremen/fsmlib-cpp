#ifndef FSM_MAIN_TESTS_H_
#define FSM_MAIN_TESTS_H_

struct TestResult {
	std::string id;
	size_t pass = 0;
	std::vector<std::string> fails;
	TestResult(const std::string &mutName) {
		id = mutName;
	}
	void printResults() {
		std::cout << "==============================================================" << std::endl;
		std::cout << "Test of " << id << std::endl;
		std::cout << "#PASS: " << pass << std::endl;
		std::cout << "#FAIL: " << fails.size() << std::endl;
		std::cout << "FAIL IDs: " << std::endl;
		for (size_t i = 0; i < fails.size(); ++i) std::cout << fails.at(i) << ", ";
		std::cout << std::endl;
		std::cout << "==============================================================" << std::endl;
	}
};

/**
 * Test Suite: Fsm::removeUnreachableNodes()
 */
void removeUnreachableNodes_TS();

/**
 * Test Suite: Fsm::transformToObservableFSM()
 */
void transformToObservableFSM_TS();

/**
 * Test Suite: Dfsm::minimise()
 */
void minimise_Dfsm_TS();

/**
 * Test Suite: Fsm::minimiseObservableFSM()
 */
void minimiseObservableFSM_TS();

/**
 * Test Suite: Fsm::minimise()
 */
void minimise_Fsm_TS();

/*
 *	Random Test Suite for test of Fsm::getCharacterisationSet().
 */
void intersect_TS();

/*
 *	Random Test Suite for test of Fsm::getCharacterisationSet().
 */
void getCharacterisationSet_Dfsm_TS();

/*
 *	Random Test Suite for test of Fsm::getCharacterisationSet().
 */
TestResult getCharacterisationSet_Fsm_TS();

/*
 *	Random Test Suite for test of FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
 *                                           const vector<shared_ptr<PkTable>>& pktblLst,
 *                                           const int maxInput)
 */
void calcDistinguishingTrace_PkTables_TS();

/*
 *	Random Test Suite for test of FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
 *                                           const vector<shared_ptr<OFSMTable>>& ofsmTblLst,
 *                                           const int maxInput,
 *                                           const int maxOutput)
 */
void calcDistinguishingTrace_OFSMTables_TS();

/*
 *	Random Test Suite for test of Fsm::calcStateIdentificationSets().
 */
void calcStateIdentificationSets_TS();

/*
 *	Random Test Suite for test of Fsm::calcStateIdentificationSetsFast().
 */
void calcStateIdentificationSetsFast_TS();

/*
 *	Random Test Suite for test of Dfsm::tMethod().
 */
void tMethod_TS();

// Test Fsm::wMethod(...)
void wMethod_Fsm_TS();

// Test Dfsm::wMethod(...)
void wMethod_Dfsm_TS();

// Test Fsm::wMethodOnMinimisedFsm(...)
void wMethodOnMinimisedFsm_TS();

// test Dfsm::wMethodOnMinimisedDfsm(...)
void wMethodOnMinimisedDfsm_TS();

// test Fsm::wpMethod(...)
void wpMethod_Fsm_TS();

// Test Dfsm::wpMethod(...)
void wpMethod_Dfsm_TS();

// Test Dfsm::wpMethodOnMinimisedDfsm(...)
void wpMethodOnMinimisedDfsm_TS();

// test Fsm::hsiMethod(...)
void hsiMethod_Fsm_TS();

// Test Dfsm::hsiMethod(...)
void hsiMethod_Dfsm_TS();

// Test Dfsm::hMethodOnMinimisedDfsm(...)
void hMethodOnMinimisedDfsm_TS();


#endif
