/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_INPUTTRACESET_H_
#define FSM_FSM_INPUTTRACESET_H_

#include <memory>
#include <unordered_set>

class InputTrace;

typedef std::unordered_set<std::shared_ptr<InputTrace>, std::hash<InputTrace>, std::equal_to<InputTrace>> InputTraceSet;

#endif //FSM_FSM_INPUTTRACESET_H_

