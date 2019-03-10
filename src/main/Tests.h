#ifndef FSM_MAIN_TESTS_H_
#define FSM_MAIN_TESTS_H_

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
void intersection_TS_Random();

/*
 *	Random Test Suite for test of Fsm::getCharacterisationSet().
 */
void getCharacterisationSet_Dfsm_TS_Random();

/*
 *	Random Test Suite for test of Fsm::getCharacterisationSet().
 */
void getCharacterisationSet_Fsm_TS_Random();

/*
 *	Random Test Suite for test of FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
 *                                           const vector<shared_ptr<PkTable>>& pktblLst,
 *                                           const int maxInput)
 */
void calcDistinguishingTrace1_TS_Random();

/*
 *	Random Test Suite for test of FsmNode::calcDistinguishingTrace(const shared_ptr<FsmNode> otherNode,
 *                                           const vector<shared_ptr<OFSMTable>>& ofsmTblLst,
 *                                           const int maxInput,
 *                                           const int maxOutput)
 */
void calcDistinguishingTrace2_TS_Random();

/*
 *	Random Test Suite for test of Fsm::calcStateIdentificationSets().
 */
void calcStateIdentificationSets_TS_Random();

/*
 *	Random Test Suite for test of Fsm::calcStateIdentificationSetsFast().
 */
void calcStateIdentificationSetsFast_TS_Random();

/*
 *	Random Test Suite for test of Dfsm::tMethod().
 */
void tMethod_TS_Random();

// Test Fsm::wMethod(...)
void wMethod_TS_Random2();

// Test Dfsm::wMethod(...)
void wMethod_Dfsm_TS_Random();

// Test Fsm::wMethodOnMinimisedFsm(...)
void wMethodOnMinimisedFsm_Fsm_TS_Random();

// test Dfsm::wMethodOnMinimisedDfsm(...)
void wMethodOnMinimisedDfsm_Dfsm_TS_Random();

// test Fsm::wpMethod(...)
void wpMethod_Fsm_TS_Random();

// Test Dfsm::wpMethod(...)
void wpMethod_Dfsm_TS_Random();

// Test Dfsm::wpMethodOnMinimisedDfsm(...)
void wpMethodOnMinimisedDfsm_Dfsm_TS_Random();

// test Fsm::hsiMethod(...)
void hsiMethod_Fsm_TS_Random();

// Test Dfsm::hsiMethod(...)
void hsiMethod_Dfsm_TS_Random();

// Test Dfsm::hMethodOnMinimisedDfsm(...)
void hMethodOnMinimisedDfsm_Dfsm_TS_Random();

#endif
