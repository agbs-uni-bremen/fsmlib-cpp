#ifndef LOGGING_H
#define LOGGING_H

#include <string>

namespace logging
{
    extern const char* fsmConversion;
    extern const char* globalLogger;
    extern const char* csvLogger;
    void initLogging();
    void setLogfileSuffix(const std::string& suffix, const std::string& loggerId = "");
}

#endif // LOGGING_H
