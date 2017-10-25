#include <string>
#include <chrono>
#include "logging/easylogging++.h"
#include "logging/Logging.h"

using namespace std;
using std::chrono::system_clock;

const char* logging::fsmConversion = "fsm-conversion";
const char* logging::globalLogger = "global-logger";

void logging::initLogging()
{
    el::Loggers::getLogger(logging::fsmConversion);
    el::Loggers::getLogger(logging::globalLogger);

    const string loggerConfigDir = "../../../src/externals/easyloggingpp_v9.95.0";
#ifdef ENABLE_DEBUG_MACRO
    //el::Configurations logConfig(loggerConfigDir + "/Debug.cfg");
    el::Loggers::configureFromGlobal((loggerConfigDir + "/Debug.cfg").c_str());
#else
    //el::Configurations logConfig(loggerConfigDir + "/Release.cfg");
    el::Loggers::configureFromGlobal((loggerConfigDir + "/Release.cfg").c_str());
#endif
    //el::Loggers::reconfigureAllLoggers(logConfig);

    const system_clock::time_point now = system_clock::now();
    std::time_t tNow = system_clock::to_time_t(now);

    char nowTextRaw[21];
    struct tm buf;
    strftime(nowTextRaw, 21, "%Y-%m-%d--%H:%M:%S", localtime_r(&tNow, &buf));
    string nowText(nowTextRaw);


    std::vector<std::string>* loggerIds = new std::vector<std::string>();
    el::Loggers::populateAllLoggerIds(loggerIds);

    for (string id : *loggerIds)
    {
        el::Logger* logger = el::Loggers::getLogger(id, false);
        string fileName = logger->configurations()->get(el::Level::Global, el::ConfigurationType::Filename)->value();
        size_t lastindex = fileName.find_last_of(".");
        if (lastindex != string::npos)
        {
            fileName = fileName.substr(0, lastindex);
        }
        fileName += "_" + nowText + ".log";
        logger->configurations()->set(el::Level::Global, el::ConfigurationType::Filename, fileName);
        logger->reconfigure();
    }
}

void logging::setLogfileSuffix(const std::string& suffix, const string& loggerId)
{
    std::vector<std::string>* loggerIds = new std::vector<std::string>();
    el::Loggers::populateAllLoggerIds(loggerIds);

    for (string id : *loggerIds)
    {
        if (id == logging::globalLogger || (loggerId != "" && id != loggerId))
        {
            continue;
        }
        el::Logger* logger = el::Loggers::getLogger(id, false);
        string fileName = logger->configurations()->get(el::Level::Global, el::ConfigurationType::Filename)->value();
        size_t lastindex = fileName.find_last_of("##");
        if (lastindex != string::npos)
        {
            fileName = fileName.substr(0, lastindex-1);
        } else {
            lastindex = fileName.find_last_of(".");
            if (lastindex != string::npos)
            {
                fileName = fileName.substr(0, lastindex);
            }
        }
        fileName += "##" + suffix + ".log";
        logger->configurations()->set(el::Level::Global, el::ConfigurationType::Filename, fileName);
        logger->reconfigure();
    }
}
