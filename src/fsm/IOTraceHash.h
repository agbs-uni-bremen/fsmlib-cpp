/*
 * Copyright. Gaël Dottel, Christoph Hilken, and Jan Peleska 2016 - 2021
 * 
 * Licensed under the EUPL V.1.1
 */
#ifndef FSM_FSM_IOTRACEHASH_H_
#define FSM_FSM_IOTRACEHASH_H_

#include <functional>

class IOTrace;

namespace std {
  template <> struct hash<IOTrace>
  {
    size_t operator()(const IOTrace& trace) const;
  };
}

#endif //FSM_FSM_IOTRACEHASH_H_

