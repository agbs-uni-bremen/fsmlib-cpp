/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_FSMLABEL_H_
#define FSM_FSM_FSMLABEL_H_

#include <memory>

#include "fsm/FsmVisitor.h"
#include "interface/FsmPresentationLayer.h"

class FsmLabel
{
private:
	/**
	The input of this label
	*/
	const int input;

	/**
	The output of this label
	*/
	const int output;

	/**
	The presentation layer used by this label
	*/
	std::shared_ptr<FsmPresentationLayer> presentationLayer;
public:
    const static int EPSILON_INPUT;
    const static int EPSILON_OUTPUT;
    const static int ERROR_OUTPUT;

	/**
	 * Create a FsmLabel
	 * @param input The input of this label
	 * @param output The ouput of this label
	 * @param maxInput The number of input
	 * @param presentationLayer The presentation layer used by this label
	 */
	FsmLabel(const int input,
             const int output,
             const std::shared_ptr<FsmPresentationLayer> presentationLayer);
    
    /**
     * Copy constructor
     */
    FsmLabel(const FsmLabel& other);

	/**
	Getter for the input
	@return The input of this label
	*/
	int getInput() const;

	/**
	Getter for the output
	@return The output of this label
	*/
	int getOutput() const;

	/**
	Check wheter or not, the 2 label are the same
	@param label1 The first label
	@param label2 The second label
	@return True if they are the same, false otherwise
	*/
	friend bool operator==(FsmLabel const & label1, FsmLabel const & label2);

    
    
    /**
     *  Accept an FsmVisitor
     */
    void accept(FsmVisitor& v);
    
    
	/**
	Check wheter or not, label1 is "smaller" than label2. this operator is needed 
	@param label1 The first label
	@param label2 The second label
	@return True if label1 is "smaller" than label2, false otherwise
	*/
	friend bool operator<(FsmLabel const & label1, FsmLabel const & label2);

	/**
	Output the FsmLabel to a standard output stream
	@param out The standard output stream to use
	@param label The FsmLabel to print
	@return The standard output stream used, to allow user to cascade <<
	*/
	friend std::ostream & operator<<(std::ostream & out, const FsmLabel & label);
};

namespace std
{
	template <>
	struct hash<FsmLabel>
	{
	public:
		size_t operator()(const FsmLabel & x) const noexcept
		{
			return ((51 + std::hash<int>()(x.getInput())) * 51 + std::hash<int>()(x.getOutput()));
		}
	};
}
#endif //FSM_FSM_FSMLABEL_H_
