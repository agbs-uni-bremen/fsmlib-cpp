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
 * Execute testsuite of Fsm::removeUnreachableNodes()
 */
TestResult removeUnreachableNodes_TS();

/**
 * Execute testsuite of Fsm::transformToObservableFSM()
 */
TestResult transformToObservableFSM_TS();

/**
 * Execute testsuite of Dfsm::minimise()
 */
TestResult minimise_Dfsm_TS();

/**
 * Execute testsuite of Fsm::minimiseObservableFSM()
 */
TestResult minimiseObservableFSM_TS();

/**
 * Execute testsuite of Fsm::minimise()
 */
TestResult minimise_Fsm_TS();

/*
 * Execute testsuite of Fsm::intersect().
 */
TestResult intersect_TS();

/*
 *	Execute testsuite of Dfsm::getCharacterisationSet().
 */
TestResult getCharacterisationSet_Dfsm_TS();

/*
 *	Execute testsuite of Fsm::getCharacterisationSet().
 */
TestResult getCharacterisationSet_Fsm_TS();

/*
 *	Execute testsuite of FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
 *                                           const vector<shared_ptr<PkTable>>& pktblLst,
 *                                           const int maxInput)
 */
TestResult calcDistinguishingTrace_PkTables_TS();

/*
 *	Execute testsuite of FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
 *                                           const vector<shared_ptr<OFSMTable>>& ofsmTblLst,
 *                                           const int maxInput,
 *                                           const int maxOutput)
 */
TestResult calcDistinguishingTrace_OFSMTables_TS();

/*
 *	Execute testsuite of Fsm::calcStateIdentificationSets().
 */
TestResult calcStateIdentificationSets_TS();

/*
 *	Execute testsuite of Fsm::calcStateIdentificationSetsFast().
 */
TestResult calcStateIdentificationSetsFast_TS();

/*
 *	Execute testsuite of Dfsm::tMethod().
 */
TestResult tMethod_TS();

// Execute testsuite of Fsm::wMethod(...)
TestResult wMethod_Fsm_TS();

// Execute testsuite of Dfsm::wMethod(...)
TestResult wMethod_Dfsm_TS();

// Execute testsuite of Fsm::wMethodOnMinimisedFsm(...)
TestResult wMethodOnMinimisedFsm_TS();

// Execute testsuite of Dfsm::wMethodOnMinimisedDfsm(...)
TestResult wMethodOnMinimisedDfsm_TS();

// Execute testsuite of Fsm::wpMethod(...)
TestResult wpMethod_Fsm_TS();

// Execute testsuite of Dfsm::wpMethod(...)
TestResult wpMethod_Dfsm_TS();

// Execute testsuite of Dfsm::wpMethodOnMinimisedDfsm(...)
TestResult wpMethodOnMinimisedDfsm_TS();

// Execute testsuite of Fsm::hsiMethod(...)
TestResult hsiMethod_Fsm_TS();

// Execute testsuite of Dfsm::hsiMethod(...)
TestResult hsiMethod_Dfsm_TS();

// Execute testsuite of Dfsm::hMethodOnMinimisedDfsm(...)
TestResult hMethodOnMinimisedDfsm_TS();


#endif
