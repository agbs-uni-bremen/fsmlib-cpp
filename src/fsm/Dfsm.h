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
#include "json/json.h"


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
	Create a DFSMTable from the DFSM
	@return The DFSMTable created
	*/
	std::shared_ptr<DFSMTable> toDFSMTable() const;
    
    std::shared_ptr<FsmPresentationLayer> createPresentationLayerFromCsvFormat(const std::string& fname);
    
    std::shared_ptr<FsmPresentationLayer> createPresentationLayerFromCsvFormat(const std::string& fname,const std::shared_ptr<FsmPresentationLayer> pl);
    
    void createDfsmTransitionGraph(const std::string& fname);
    
public:
	/**
	Create a DFSM from a file description
	@param fname The name of the file containing the FSM informations
	@param fsmName The name of the DFSM
	@param maxNodes Number of DFSM states
	@param maxInput Maximal value of the (integer) input alphabet - admissible
	values are 0..maxInput
	@param maxOutput Maximal value of (integer) output alphabet - admissible
	values are 0..maxOutput
	@param presentationLayer The presentation layer used by the DFSM
     
     \note This constructor is deprecated and will be removed in a future version
	*/
	Dfsm(const std::string & fname,
         const std::string & fsmName,
         const int maxNodes,
         const int maxInput,
         const int maxOutput,
         const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    
    /**
     Create a DFSM from a file description
     @param fname The name of the file containing the FSM informations
     @param presentationLayer The presentation layer used by the DFSM
     @param fsmName The name of the DFSM
     
     The parameters maxNodes, maxInput, and maxOutput used in the constructor above
     are determined from the input file specifying the FSM.
     */
    Dfsm(const std::string & fname,
         const std::shared_ptr<FsmPresentationLayer> presentationLayer,
         const std::string & fsmName);

    
    /**
     Create a random mutant of this DFSM
     */
    void createAtRandom();
    
	/**
	Random creation of a completely defined deterministic FSM
	@param fsmName The name of the DFSM
	@param maxNodes Number of DFSM states
	@param maxInput Maximal value of the (integer) input alphabet - admissible
	values are 0..maxInput
	@param maxOutput Maximal value of (integer) output alphabet - admissible
	values are 0..maxOutput
	@param presentationLayer The presentation layer used by the DFSM
	*/
	Dfsm(const std::string & fsmName,
         const int maxNodes,
         const int maxInput,
         const int maxOutput,
         const std::shared_ptr<FsmPresentationLayer> presentationLayer);

	/**
	Create a DFSM from a list of nodes
	@param fsmName The name of the DFSM
	@param maxInput Maximal value of the (integer) input alphabet - admissible
	values are 0..maxInput
	@param maxOutput Maximal value of (integer) output alphabet - admissible
	values are 0..maxOutput
	@param lst The list of nodes of the DFSM
	@param presentationLayer The presentation layer used by the DFSM
	*/
	Dfsm(const std::string & fsmName,
         const int maxInput,
         const int maxOutput,
         const std::vector<std::shared_ptr<FsmNode>> lst,
         const std::shared_ptr<FsmPresentationLayer> presentationLayer);

	/**
	Create a DFSM from the equivalent deterministic FSM
	@param fsm The deterministic FSM
	*/
	Dfsm(const Fsm & fsm);
    
    
    /**
     * Create a DFSM from a file written in csv format.
     * @item fields are separated by ;
     * @item First row starts with empty field, followed by
     *       input event names to be used by presentation layer
     * @item Consecutive rows are ordered by DFSM states,
     *       starting with the initial state
     * @item Each of these rows contain the state name in the first 
     *       column, as it should be used by the presentation layer.
     * @item In row indexed by state s, column indexed by input x
     *       contains entry 
     *                  s' / y
     *       where s' is the post-state reached when applying x in state s,
     *       and y is the associated output.
     *
     * @param fname Filename with DFSM specification in csv format.
     * @param fsmName Name of the DFSM to be used in output files etc.
     */
    Dfsm(const std::string& fname,
         const std::string& fsmName);
    
    /**
     *  This constructor acts like the previous, but uses a
     *  a given presentation layer. Typically, this is required
     *  if an implementation DFSM is read from a file where
     *  inputs and outputs occur in a different order, and
     *  where some inputs/outputs occurring in the reference model
     *  do not occur in the implementation. Also, the implementation
     *  might use additional inputs and/or outputs that do not
     *  occur in the reference model. These are appended to the respective
     *  lists of the given presentation layer.
     *
     *  The constructor ensures that both implementation model and
     *  reference model use the same internal numbers for the
     *  same inputs and outputs.
     */
    Dfsm(const std::string& fname,
         const std::string& fsmName,
         const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    
    
    /**
     *  Construct an DFSM from a json model file. The DFSM is completely
     *  specified: in every state of the json model, inputs x that do
     *  not occur in any outgoing transition, are associated with
     *  self-loop transitions that have an additional NO OPERATION
     *  action as output. The NOP action is added as last action
     *  to the output alphabet.
     */
    Dfsm(const Json::Value& jsonModel);
    
    /**
     *  This constructor acts like the previous, but uses a 
     *  a given presentation layer. Typically, this is required
     *  if an implementation DFSM is read from a file where
     *  inputs and outputs occur in a different order, and 
     *  where some inputs/outputs occurring in the reference model
     *  do not occur in the implementation. Also, the implementation
     *  might use additional inputs and/or outputs that do not
     *  occur in the reference model. These are appended to the respective 
     *  lists of the given presentation layer.
     *
     *  The constructor ensures that both implementation model and
     *  reference model use the same internal numbers for the
     *  same inputs and outputs.
     */
    Dfsm(const Json::Value& jsonModel,
         const std::shared_ptr<FsmPresentationLayer> presentationLayer);


	/**
	Minimise this DFSM
	As a side effect, create the DFSM table and all Pk tables needed.
	@return the new minimised DFSM
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
	@return Characterisation set
	*/
	IOListContainer getCharacterisationSet();

	/**
	Apply input trace to the initial state of the FSM.
	\note This operation is only applicable to deterministic FSMs
	@param i input trace
	@return IOTrace containing the maximal sub-trace of i that could be processed as inputs
	(if FSM is completely defined this is always the complete trace i) and the associated
	outputs of the same length (the outputs are uniquely determined by the inputs, since the FSM
	is deterministic.
	*/
	IOTrace applyDet(const InputTrace & i);

	/**
	Check whether IOTrace is in the language of the DFSM
	@param io IOTrace to be checked against the DFSM
	@return True if io is in the language of the DFSM,
	False otherwise
	*/
	bool pass(const IOTrace & io);

   /**
	* Perform test generation by means of the W-Method.
    * The DFSM this method is applied to is regarded as the reference
    * model. If an implementation DFSM passes this test suite, this proves
    * language equivalence (I/O-equivalence) between reference DFSM and
    * implementation DFSM, provided that the implementation DFSM in minimised
    * form does not have more than numAddStates additional states, when compared
    * to the minimised reference DFSM. If this assumption does not hold,
    * the test suite may not uncover certain errors in the implementation DFSM.
    *
    * The reference DFSM will first be minimised, and then the
    * proper W-Method is applied to the minimised DFSM.
	* @param numAddStates The maximal number of additional states,
    *                     which the implementation DFSM in minimised 
    *                     for may have, when compared to the reference
    *                     model in minimised form.
	* @return A test suite
    *
    * @note The size of the test suite to be produced grows exponentially
    *        with numAddStates
    *
    * @note If it is already known that the reference DFSM is minimal,
    *       then method wMethodOnMinimisedDfsm() should rather be used,
    *       since it avoids unnecessary minisation steps.
	*/
	IOListContainer wMethod(const unsigned int numAddStates);

    /**
     *  Apply the W-Method on a DFSM that is already minimised
     */
    IOListContainer wMethodOnMinimisedDfsm(const unsigned int numAddStates);


	/**
	* Perform test generation by means of the Wp Method. The algorithm
	* we have implemented is applicable to both nondeterministic and
	* deterministic FSMs. It relies, however, on the existence of
	* OFSM tables which we only calculate for nondeterministic FSMs.
	* Therefore we need a wrapper method for DFSMs which first
	* calculates OFSM tables (by means of a call to minimiseObervableFSM())
	* and then calls the wpMethod() operation of the super class Fsm.
     * @param numAddStates The maximal number of additional states,
     *                     which the implementation DFSM in minimised
     *                     form may have, when compared to the reference
     *                     model in minimised form.
	* @return A test suite
	*/
	IOListContainer wpMethod(const unsigned int numAddStates);

    /**
     * WORK IN PROGRESS
     * Perform test generation by means of the HSI-Method. The algorithm
     * we have implemented is applicable to both nondeterministic and
     * deterministic FSMs.
     * @param numAddStates The maximal number of additional states,
     *                     which the implementation DFSM in minimised
     *                     form may have, when compared to the reference
     *                     model in minimised form.
     * @return a test suite
     */
    IOListContainer hsiMethod(const unsigned int numAddStates);
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
    
    /**
     *  Output DFSM in tabular format as *.csv file
     *  The table consists of one row per state and one column per input.
     *  The table entry for state s and input x consists of the pair
     *  s'/y, where s' is the post-state of the transition triggered
     *  by x in state s, and y is the corresponding output.
     *  @param fname Name of the output file without extension; 
     *         extension .csv is added internally to the filename.
     */
    void toCsv(const std::string& fname);
    
};
#endif //FSM_FSM_DFSM_H_
