/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_FSMTRANSITION_H_
#define FSM_FSM_FSMTRANSITION_H_

#include <memory>

#include "fsm/FsmLabel.h"

class FsmNode;

class FsmTransition
{
private:
	/**
	 *The node from which the transition comes
	 */
	std::shared_ptr<FsmNode> source;

	/**
	 * The node where the transition go
	 */
	std::shared_ptr<FsmNode> target;

	/**
	 * The label of this transition
	 */
	std::shared_ptr<FsmLabel> label;
public:
	/**
	Create a FsmTransition
	\param source The node from which the transition come
	\param target The node where the transition go
	\param label The label of this transition
	*/
	FsmTransition(const std::shared_ptr<FsmNode>  source,
                  const std::shared_ptr<FsmNode>  target,
                  const std::shared_ptr<FsmLabel> label);

	/**
	 * Getter for the source
	 * @return The source node of the transition
	 */
	std::shared_ptr<FsmNode> getSource();

	/**
	Getter for the target
	@return The target node of tghis transition
	*/
	std::shared_ptr<FsmNode> getTarget();
    
    /** setter functions */
    void setSource(std::shared_ptr<FsmNode> src);
    void setTarget(std::shared_ptr<FsmNode> tgt);
    void setLabel(std::shared_ptr<FsmLabel> lbl);

    
	/**
	Getter for the label
	@return The label of this transition
	*/
    std::shared_ptr<FsmLabel> getLabel();

	/**
	Output the FsmTransition to a standard output stream
	\param out The standard output stream to use
	\param transition The FsmTransition to print
	@return The standard output stream used, to allow user to cascade <<
	*/
	friend std::ostream & operator<<(std::ostream & out, FsmTransition & transition);
};
#endif //FSM_FSM_FSMTRANSITION_H_
