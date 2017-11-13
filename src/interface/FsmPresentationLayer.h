
/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_INTERFACE_FSMPRESENTATIONLAYER_H_
#define FSM_INTERFACE_FSMPRESENTATIONLAYER_H_

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

class FsmPresentationLayer
{
private:
	/**
	 * A vector containing a string for each input
	 */
	std::vector<std::string> in2String;

	/**
	 * A vector containing a string for each output
	 */
	std::vector<std::string> out2String;

	/**
	 * A vector containing a string for each state
	 */
	std::vector<std::string> state2String;
public:
    const static std::string FAIL_STATE_NAME;
    const static std::string FAIL_OUTPUT;
	/**
	 * Create a new presentation layer containing nothing
	 */
	FsmPresentationLayer();
    
    /**
     * Copy constructor
     */
    FsmPresentationLayer(const FsmPresentationLayer& pl);

	/**
	 * Create a new presentation layer
	 * @param in2String A vector containing a string for each input
	 * @param out2String A vector containing a string for each output
	 * @param state2String A vector containing a string for each state
	 */
	FsmPresentationLayer(const std::vector<std::string>& in2String,
                         const std::vector<std::string>& out2String,
                         const std::vector<std::string>& state2String);

	/**
	 * Create a new presentation layer
	 * @param inputs A file name, in which each line stand for one input
	 * @param outputs A file name, in which each line stand for one output
	 * @param states A file name, in which each line stand for one state
	 */
	FsmPresentationLayer(const std::string & inputs,
                         const std::string & outputs,
                         const std::string & states);

    static std::shared_ptr<FsmPresentationLayer> mergeAlphabets(std::shared_ptr<const FsmPresentationLayer> plA, std::shared_ptr<const FsmPresentationLayer> plB);

    void setState2String(std::vector<std::string> state2String);
    void addState2String(std::string name);
    void removeState2String(const int index);
    int addOut2String(std::string name);
    int addOut2String(const int i, std::string name);
    int addIn2String(std::string name);
    int addIn2String(const int i, std::string name);
    void truncateState2String(const int index);
    void truncateIn2String(const int index);
    void truncateOut2String(const int index);

	/**
	 * Getter for a particular input name
	 * @param id The id of the input
	 * @return The name if this input
	 */
	std::string getInId(const unsigned int id) const;

	/**
	 * Getter for a particular output name
	 * @param id The id of the output
	 * @return The name if this output
	 */
	std::string getOutId(const unsigned int id) const;

	/**
	 * Getter for a particular state name
	 * @param id The id of the state
	 * @param prefix If the prefix is not empty, it will be added at he beginning of the name before returning it
	 * @return The name if this state
	 */
	std::string getStateId(const unsigned int id, const std::string & prefix) const;
    
    /**
     *  Get a copy of the in2string vector
     */
    std::vector<std::string> getIn2String() { return in2String; }
    
    /**
     *  Get a copy of the out2string vector
     */
    std::vector<std::string> getOut2String() { return out2String; }
    
    /**
     *  Get a copy of the state2string vector
     */
    std::vector<std::string> getState2String() { return state2String; }
    
    /**
     *  Convert input name to input number
     */
     int in2Num(const std::string& name);
    
    /**
     *  Convert output name to output number
     */
     int out2Num(const std::string& name);
    
    /**
     *  Convert state name to state number
     */
     int state2Num(const std::string& name);

	/**
	 * Dump the current inputs into an output stream
	 * @param out The standard output to use
	 */
	void dumpIn(std::ostream & out) const;

	/**
	 * Dump the current outputs into an output stream
	 * @param out The standard output to use
	 */
	void dumpOut(std::ostream & out) const;

	/**
	 * Dump the current states into an output stream
	 * @param out The standard output to use
	 */
	void dumpState(std::ostream & out) const;

	/**
	 * Compare two presentation layer to check if they are the same or not
	 * @param otherPresentationLayer The other presentation layer to be compared
	 */
	bool compare(std::shared_ptr<FsmPresentationLayer> otherPresentationLayer);
    FsmPresentationLayer& operator=(FsmPresentationLayer& other);
};
#endif //FSM_INTERFACE_FSMPRESENTATIONLAYER_H_
