#ifndef LOGGING_H
#define LOGGING_H

#include <string>

namespace logging
{
    extern const char* fsmConversion;
    void initLogging();
    void setLogfileSuffix(const std::string& suffix);
}

#endif // LOGGING_H
