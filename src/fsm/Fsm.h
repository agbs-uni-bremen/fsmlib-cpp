/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_FSM_H_
#define FSM_FSM_FSM_H_

#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

class Dfsm;
class FsmNode;
class Tree;
class OutputTree;
class InputTrace;
class FsmPresentationLayer;
class OFSMTable;
class IOListContainer;
class TestSuite;

enum Minimal
{
	True, False, Maybe
};

class Fsm
{
protected:
    
    /**
     *  default constructor without effect - needed by sub-classes
     */
    Fsm(const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    


    
	std::string name;
	std::vector<std::shared_ptr<FsmNode>> nodes;
	std::shared_ptr<FsmNode> currentParsedNode;
	int maxInput;
	int maxOutput;
    int maxState;
	int initStateIdx;
	std::vector<std::shared_ptr<OFSMTable>> ofsmTableLst;
	std::shared_ptr<Tree> characterisationSet;
	std::vector<std::shared_ptr<Tree>> stateIdentificationSets;
	std::shared_ptr<FsmPresentationLayer> presentationLayer;
	Minimal minimal;
	std::shared_ptr<FsmNode> newNode(const int id, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>> p);
	bool contains(const std::vector<std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>>>& lst, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>> p);
	bool contains(const std::vector<std::shared_ptr<FsmNode>>& lst, const std::shared_ptr<FsmNode> n);
	std::shared_ptr<FsmNode> findp(const std::vector<std::shared_ptr<FsmNode>>& lst, const std::shared_ptr<std::pair<std::shared_ptr<FsmNode>, std::shared_ptr<FsmNode>>> p);
	void parseLine(const std::string & line);
	void readFsm(const std::string & fname);
    
    void parseLineInitial (const std::string & line);
    void readFsmInitial (const std::string & fname);

    
	std::string labelString(std::unordered_set<std::shared_ptr<FsmNode>>& lbl) const;
public:

    /**
     *  Constructor creating an FSM from file - used only internally
     */
    Fsm(const std::string & fname,
        const std::string & fsmName,
        const int maxNodes,
        const int maxInput,
        const int maxOutput,
        const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    
    /** 
     *  Constructor creating an FSM specified in a file.
     *  \param fname Filename of a text file (typically with extension *.fsm),
     *  where each row is formatted as
     *       pre-state input output post-state
     * \item pre-state is a number in range 0..(number of states -1) specifying
     *       the FSM pre-state
     * \item input is a number in range 0..maxInput specifying the input applied to
     *       the pre-state
     * \item output is a number in range 0..maxOutput specifying the output
     *       produced by the FSM when transiting to post-state on 'input'
     * \item post-state is a number in range 0..(number of states -1)
     * The initial state is expected to be the pre-state of the first row.
     * The parameters 'number of states', maxInput, and maxOutput are
     * determined when parsing the input file fname.
     *
     * \param presentationLayer Pointer to instance of a presentation layer
     *        which associates each input number in range 0..maxInput with
     *        an input name, each output number in range 0..maxOutput with
     *        an output name, and each state number in range 0..(number of states -1)
     *        with a state name.
     *
     * \param fsmName name of the FSM
     *
     */
    Fsm(const std::string& fname,
        const std::shared_ptr<FsmPresentationLayer> presentationLayer,
        const std::string& fsmName);
    
    
    
    
    

	/**
	Constructor for creating an FSM from a list of FsmNodes that have
     been created beforehand, together with the outgoing transitions of each FsmNode/.
	\param fsmName  name of the FSM (appears in every node)
	\param maxInput maximal value of the (integer) input alphabet - admissible
	values are 0..maxInput
	\param maxOutput maximal value of (integer) output alphabet - admissible
	values are 0..maxOutput
	*/
	Fsm(const std::string & fsmName,
        const int maxInput,
        const int maxOutput,
        const std::vector<std::shared_ptr<FsmNode>> lst,
        const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    
    
    
    
    /**
     * Write FSM to text file using the standard format in each line
     * which is also used for instantiating an FSM from file:
     *    pre-state input output post-state
     * \item pre-state is a number in range 0..(number of states -1) specifying
     *       the FSM pre-state
     * \item input is a number in range 0..maxInput specifying the input applied to
     *       the pre-state
     * \item output is a number in range 0..maxOutput specifying the output
     *       produced by the FSM when transiting to post-state on 'input'
     * \item post-state is a number in range 0..(number of states -1)
     */
	void dumpFsm(std::ofstream & outputFile) const;
	std::shared_ptr<FsmNode> getInitialState() const;

	/**
	This is the getter for the name of the FSM
	\return Name of the FSM
	*/
	std::string getName() const;
	virtual int getMaxNodes() const;//TODO NOT PRESENT IN JAVA
	int getMaxInput() const;//TODO NOT PRESENT IN JAVA
	int getMaxOutput() const;//TODO NOT PRESENT IN JAVA
	std::vector<std::shared_ptr<FsmNode>> getNodes() const;//TODO NOT PRESENT IN JAVA
	std::shared_ptr<FsmPresentationLayer> getPresentationLayer() const;//TODO NOT PRESENT IN JAVA
	int getInitStateIdx() const;//TODO NOT PRESENT IN JAVA
	void resetColor();
	void toDot(const std::string & fname);

	/**
	Create a new FSM that represents the intersection of this and the other FSM
	\param f the other FSM
	\return a new FSM which equals the intersection of this and f
	*/
	Fsm intersect(const Fsm & f);
    
    /**
     * Generate the state cover of an arbitrary FSM
     */
	std::shared_ptr<Tree> getStateCover();
    
    /**
     * Generate the transition cover of an arbitrary FSM
     */
	std::shared_ptr<Tree> getTransitionCover();
    
    /**
     *  Apply an input trace to an FSM and return its
     *  resulting output tree.
     */
	OutputTree apply(const InputTrace & itrc);
    
    /**
     *  Transform an FSM to its observable equivalent.
     */
	Fsm transformToObservableFSM() const;

	/**
	Check this FSM with respect to observability
	\return true if and only if the FSM is observable
	*/
	bool isObservable() const;
	Minimal isMinimal() const;

	/**
	Create the minimal observable FSM which is equivalent to this FSM.
	\pre This method can only be applied to an observable OFSM
	\return minimal observable FSM which is equivalent to this FSM
	*/
	Fsm minimiseObservableFSM();

	/**
	Create the minimal observable FSM which is equivalent to this FSM.
	If this FSM is not observable, an observable equivalent is created
	first. The observable FSM is then minimised.
	\return minimal observable FSM which is equivalent to this FSM
	*/
	Fsm minimise();
	bool isCharSet(const std::shared_ptr<Tree> w) const;
	void minimiseCharSet(const std::shared_ptr<Tree> w);

	/**
	Calculate the characterisation set W of a (possibly nondeterministic) FSM.
	In addition, calculate the state identification
	sets W_i containing prefixes of input traces
	of W, such that W_i distinguishes FSM state number i from every other
	FSM state.
	\pre The FSM must be observable and minimised
	\return The characterisation set W, represented by an instance of
	IOListContainer
	\note As a side effect, the state identification sets W_i
	are stored in the list stateIdentificationSets.
	*/
	IOListContainer getCaracterisationSet();

	/**
	Calculate the state identification sets. The sets are stored
	in list stateIdentificationSets, ordered by the FSM state numbers.
	\pre The FSM must be observable and minimal, and the characterisation
	set must have been previously calculated using operation
	getCharacterisationSet().
	*/
	void calcStateIdentificationSets();
	void appendStateIdentificationSets(const std::shared_ptr<Tree> Wp2) const;

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
	TestSuite createTestSuite(const IOListContainer & testCases);
	bool isCompletelyDefined() const;

	/**
	Check if FSM is deterministic
	\return true if FSM is deterministic
	*/
	bool isDeterministic() const;
	void setPresentationLayer(const std::shared_ptr<FsmPresentationLayer> ppresentationLayer);//TODO NOT PRESENT IN JAVA
	friend std::ostream & operator<<(std::ostream & out, const Fsm & fsm);
};
#endif //FSM_FSM_FSM_H_
