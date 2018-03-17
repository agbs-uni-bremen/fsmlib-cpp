#include <string>
#include <chrono>
#include "logging/easylogging++.h"
#include "logging/Logging.h"

using namespace std;
using std::chrono::system_clock;

const char* logging::fsmConversion = "fsm-conversion";
const char* logging::globalLogger = "global-logger";
const char* logging::csvLogger = "csv-logger";

void logging::initLogging()
{
    el::Loggers::getLogger(logging::fsmConversion);
    el::Loggers::getLogger(logging::globalLogger);
    el::Loggers::getLogger(logging::csvLogger);

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
    strftime(nowTextRaw, 21, "%Y-%m-%d--%H-%M-%S", localtime_r(&tNow, &buf));
    string nowText(nowTextRaw);


    std::vector<std::string>* loggerIds = new std::vector<std::string>();
    el::Loggers::populateAllLoggerIds(loggerIds);

    for (string id : *loggerIds)
    {
        el::Logger* logger = el::Loggers::getLogger(id, false);
        string fileName = logger->configurations()->get(el::Level::Global, el::ConfigurationType::Filename)->value();
        string fileExtension = "";
        size_t lastindex = fileName.find_last_of(".");
        if (lastindex != string::npos)
        {
            fileExtension = fileName.substr(lastindex);
            fileName = fileName.substr(0, lastindex);
        }
        fileName += "_" + nowText + fileExtension;
        logger->configurations()->set(el::Level::Global, el::ConfigurationType::Filename, fileName);
        logger->reconfigure();
    }
}

void logging::setLogfileSuffix(const std::string& suffix, const string& loggerId)
{
    std::vector<std::string>* loggerIds = new std::vector<std::string>();
    el::Loggers::populateAllLoggerIds(loggerIds);
    el::base::LogStreamsReferenceMap* streams = el::base::elStorage->registeredLoggers()->logStreamsReference();

    for (string id : *loggerIds)
    {
        if ((id == logging::globalLogger && loggerId != logging::globalLogger)
                || (id == logging::csvLogger && loggerId != logging::csvLogger)
                || (loggerId != "" && id != loggerId))
        {
            continue;
        }
        el::Logger* logger = el::Loggers::getLogger(id, false);
        string fileName = logger->configurations()->get(el::Level::Global, el::ConfigurationType::Filename)->value();
        string fileExtension = "";
        size_t poundIndex = fileName.find_last_of("##");
        size_t dotIndex = fileName.find_last_of(".");
        if (dotIndex != string::npos)
        {
            fileExtension = fileName.substr(dotIndex);
        }
        if (poundIndex != string::npos)
        {
            fileName = fileName.substr(0, poundIndex-1);
        } else {
            size_t dotIndex = fileName.find_last_of(".");
            if (dotIndex != string::npos)
            {
                fileName = fileName.substr(0, dotIndex);
            }
        }
        fileName += "##" + suffix + fileExtension;

        el::Configurations newConf;
        newConf.set(el::Level::Global, el::ConfigurationType::Filename, fileName);
        logger->configure(newConf);

        for (auto it = streams->begin(); it != streams->end(); ++it)
        {
            if (it->first == fileName)
            {
                streams->erase(it);
                break;
            }
        }
    }
}
