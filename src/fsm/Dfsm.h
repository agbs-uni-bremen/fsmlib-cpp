/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_DFSM_H_
#define FSM_FSM_DFSM_H_

#include <stdlib.h>
#include <time.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "fsm/DFSMTable.h"
#include "fsm/Fsm.h"

class PkTable;
class IOTrace;

class Dfsm : public Fsm
{
private:
	//TODO
	std::shared_ptr<DFSMTable> dfsmTable;

	//TODO
	std::vector<std::shared_ptr<PkTable>> pktblLst;

	/**
	Create a random DFSM
	*/
	void createAtRandom();

	/**
	Create a DFSMTable from the DFSM
	\return The DFSMTable created
	*/
	std::shared_ptr<DFSMTable> toDFSMTable() const;
    
    
public:
	/**
	Create a DFSM from a file description
	\param fname The name of the file containing the FSM informations
	\param fsmName The name of the DFSM
	\param maxNodes Number of DFSM states
	\param maxInput Maximal value of the (integer) input alphabet - admissible
	values are 0..maxInput
	\param maxOutput Maximal value of (integer) output alphabet - admissible
	values are 0..maxOutput
	\param presentationLayer The presentation layer used by the DFSM
     
     \note This constructor is deprecated and will be removed in a furture version
	*/
	Dfsm(const std::string & fname,
         const std::string & fsmName,
         const int maxNodes,
         const int maxInput,
         const int maxOutput,
         const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    
    /**
     Create a DFSM from a file description
     \param fname The name of the file containing the FSM informations
     \param presentationLayer The presentation layer used by the DFSM
     \param fsmName The name of the DFSM
     
     The parameters maxNodes, maxInput, and maxOutput used in the constructor above
     are determined from the input file specifying the FSM.
     */
    Dfsm(const std::string & fname,
         const std::shared_ptr<FsmPresentationLayer> presentationLayer,
         const std::string & fsmName);

	/**
	Random creation of a completely defined deterministic FSM
	\param fsmName The name of the DFSM
	\param maxNodes Number of DFSM states
	\param maxInput Maximal value of the (integer) input alphabet - admissible
	values are 0..maxInput
	\param maxOutput Maximal value of (integer) output alphabet - admissible
	values are 0..maxOutput
	\param presentationLayer The presentation layer used by the DFSM
	*/
	Dfsm(const std::string & fsmName,
         const int maxNodes,
         const int maxInput,
         const int maxOutput,
         const std::shared_ptr<FsmPresentationLayer> presentationLayer);

	/**
	Create a DFSM from a list of nodes
	\param fsmName The name of the DFSM
	\param maxInput Maximal value of the (integer) input alphabet - admissible
	values are 0..maxInput
	\param maxOutput Maximal value of (integer) output alphabet - admissible
	values are 0..maxOutput
	\param lst The list of nodes of the DFSM
	\param presentationLayer The presentation layer used by the DFSM
	*/
	Dfsm(const std::string & fsmName,
         const int maxInput,
         const int maxOutput,
         const std::vector<std::shared_ptr<FsmNode>> lst,
         const std::shared_ptr<FsmPresentationLayer> presentationLayer);

	/**
	Create a DFSM from the equivalent deterministic FSM
	\param fsm The deterministic FSM
	*/
	Dfsm(const Fsm & fsm);

	/**
	Minimise this DFSM
	As a side effect, create the DFSM table and all Pk tables needed.
	\return the new minimised DFSM
	*/
	Dfsm minimise();

	/**
	Output the DFSM table and the Pk-tables in LaTeX format.
	The output file name is hard-coded as tables.tex.
	\note This requires that the operation minimise() has been called before.
	*/
	void printTables() const;

	/**
	Calculate characterisation set according to Gill's algorithm
	based on Pk-tables
	\pre the DFSM has already been minimised, and it is completely defined
	\return Characterisation set
	*/
	IOListContainer getCharacterisationSet();

	/**
	Apply input trace to the initial state of the FSM.
	\note This operation is only applicable to deterministic FSMs
	\param i input trace
	\return IOTrace containing the maximal sub-trace of i that could be processed as inputs
	(if FSM is completely defined this is always the complete trace i) and the associated
	outputs of the same length (the outputs are uniquely determined by the inputs, since the FSM
	is deterministic.
	*/
	IOTrace applyDet(const InputTrace & i);

	/**
	Check whether IOTrace is in the language of the DFSM
	\param io IOTrace to be checked against the DFSM
	\return True if io is in the language of the DFSM,
	False otherwise
	*/
	bool pass(const IOTrace & io);

	/**
	Perform test generation by means of the W Method.
	\param m Maximum number of states
	\return A test suite
	*/
	IOListContainer wMethod(const unsigned int m);


	/**
	Perform test generation by means of the Wp Method. The algorithm
	we have implemented is applicable to both nondeterministic and
	deterministic FSMs. It relies, however, on the existence of
	OFSM tables which we only calculate for nondeterministic FSMs.
	Therefore we need a wrapper method for DFSMs which first
	calculates OFSM tables (by means of a call to minimiseObervableFSM())
	and then calls the wpMethod() operation of the super class Fsm.
	\param m Maximum number of states
	\return A test suite
	*/
	IOListContainer wpMethod(const int m);
    
    /**
     * Perform test generation by means of the T-Method. The algorithm 
     * applies to deterministic FSMs which are completely defined.
     * It produces test cases from the transition cover. The resulting 
     * test suite is complete for I/O-equivalence in the fault domain
     * of all completely defined DFSMs M' differing from the reference
     * model M only by zero or more output faults. The number of in M and M'
     * must be identical, and the transition arrows and input labels
     * must be identical as well to ensure completeness.
     */
    IOListContainer tMethod();
    
};
#endif //FSM_FSM_DFSM_H_
