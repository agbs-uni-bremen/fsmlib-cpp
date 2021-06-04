#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <unordered_set>
#include <map>
#include <fstream>
#include <memory>
#include <vector>
#include <stack>



typedef enum {
    EQUIVALENCE,
    REDUCTION
} conformance_relation_t;

static size_t k; // number of times to repeat application of each test
static conformance_relation_t conformanceRelation;
static std::string testSuiteFile;


/** functions expected from the wrapper */
extern void sut_init();
extern void sut_reset();
extern const std::string sut(const std::string input);



/** definitions for a simple type of tree to represent the test suite with */
typedef std::pair<std::string,std::string> IOPair;

struct TestNode {
    bool visited = false;
    std::shared_ptr<TestNode> parent = nullptr;
    IOPair pairFromParent;
    std::map<IOPair,std::shared_ptr<TestNode>> nexts {};
};



/**
 * Write program usage to standard error.
 * @param name program name as specified in argv[0]
 */
static void printUsage(char* name) {
    std::cerr << "usage: " 
              << name
              << "[-k input application repetitions] [-equivalence|-reduction] test-suite-file " << std::endl;
}




static void parseParameters(int argc, char* argv[]) {
    
    // set default values
    k = 100;
    conformanceRelation = EQUIVALENCE;

    bool testSuiteFileNameExists = false;    
    
    for ( int p = 1; p < argc; p++ ) {
        
        if ( strcmp(argv[p],"-k") == 0 ) {
            if ( argc < p+2 ) {
                std::cerr << argv[0] << ": missing number of input application repetitions" << std::endl;
                printUsage(argv[0]);
                exit(1);
            }
            else {
                k = atoi(argv[++p]);
            }
        }
        else if ( strcmp(argv[p],"-reduction") == 0 ) {
            conformanceRelation = REDUCTION;
        }
        else if ( strcmp(argv[p],"-equivalence") == 0 ) {
            conformanceRelation = EQUIVALENCE;
        }
        else {
            if (testSuiteFileNameExists) {
                std::cerr << argv[0] << ": illegal parameter `" << argv[p] << "'" << std::endl;
                printUsage(argv[0]);
                exit(1);
            } else {
                testSuiteFile = argv[p];
                testSuiteFileNameExists = true;
            }
        }        
    }
    
    if ( !testSuiteFileNameExists ) {
        std::cerr << argv[0] << ": missing test suite file" << std::endl;
        printUsage(argv[0]);
        exit(1);
    }    
}

static std::vector<IOPair> parseTestCase(std::string& tc) {
    std::vector<IOPair> result;
    std::string s = tc;
    size_t pos = s.find("(");
    while (pos != std::string::npos) {        
        s.erase(0,pos+1);
        pos = s.find("/");
        if (pos == std::string::npos) {
            std::cerr << "incorrect test case format: " << tc << std::endl;
            exit(1);
        }
        std::string x = s.substr(0,pos);
        s.erase(0,pos+1);
        pos = s.find(")");
        if (pos == std::string::npos) {
            std::cerr << "incorrect test case format: " << tc << std::endl;
            exit(1);
        }
        std::string y = s.substr(0,pos);
        s.erase(0,pos+1);
        pos = s.find("(");

        IOPair xy = std::make_pair(x,y);
        result.push_back(xy);
    }
    return result;
}

static void addTestCase(std::shared_ptr<TestNode> root, std::string& tc) {

    std::shared_ptr<TestNode> node = root;
    std::vector<IOPair> parsedTC = parseTestCase(tc);
    
    for (auto& xy : parsedTC) {
        if (node->nexts.count(xy) == 0) {
            auto target = std::make_shared<TestNode>();
            target->parent = node;
            target->pairFromParent = xy;
            node->nexts[xy] = target;
        }
        node = node->nexts[xy];
    }
}

static void printTestTree(std::shared_ptr<TestNode> root, size_t depth = 1) {
    std::cout << (root->visited ? "V" : "N");

    if (root->nexts.empty()) {
        std::cout << std::endl;
    } else {
        bool isFirst = true;
        std::cout<<"\t";
        for (auto& entry : root->nexts) {
            if (!isFirst) {
                std::cout << std::string(depth,'\t');
            }
            isFirst = false;
            std::cout << entry.first.first << "/" << entry.first.second;
            printTestTree(entry.second, depth+1);
        }
    }
}

static std::shared_ptr<TestNode> generateTestSuiteTree() {
    std::shared_ptr<TestNode> testSuite = std::make_shared<TestNode>();

    std::ifstream file(testSuiteFile);
    std::string testCaseLine;
    while (std::getline(file, testCaseLine))
    {
        addTestCase(testSuite,testCaseLine);
    }
    return testSuite;
}

static std::stack<IOPair> getPathToNode(std::shared_ptr<TestNode> node) {
    std::stack<IOPair> result;
    std::shared_ptr<TestNode> currentNode = node;
    if (currentNode == nullptr) 
        return result;

    while (currentNode->parent != nullptr) {
        result.push(currentNode->pairFromParent);
        currentNode = currentNode->parent; 
    }
    return result;
}

static void printPathToNode(std::shared_ptr<TestNode> node) {
    std::stack<IOPair> st = getPathToNode(node);
    bool isFirst = true;
    while (!st.empty()) {
        auto& xy = st.top();
        st.pop();
        if (!isFirst) std::cout << ".";
        std::cout << "(" << xy.first << "/" << xy.second << ")";
        isFirst = false;
    }
}

static void applyTestCase(std::string tcString, std::vector<IOPair>& tc,std::shared_ptr<TestNode>testSuiteTree) {
    std::shared_ptr<TestNode> node = testSuiteTree;

    sut_reset();

    for (auto& xy : tc) {
        std::string response = sut(xy.first);

        IOPair sutPair = std::make_pair(xy.first, response);
        if (node->nexts.count(sutPair) == 0) {
            std::cerr << "Failure observed for test case " << tcString << std::endl;
            std::cerr << "  observed invalid response " << response << " to input " << xy.first << " after prefix ";
            printPathToNode(node);
            std::cerr << std::endl;
            std::cerr << "SUT is does not pass the test suite." << std::endl;
            exit(1);
        }

        node = node->nexts[sutPair];
        node->visited = true;

        // continue only if the observed response matches the response for this test case
        if (sutPair.second != xy.second) return;        
    }
}

static void applyTestSuite(std::shared_ptr<TestNode> testSuiteTree) {
    std::ifstream file(testSuiteFile);
    std::string testCaseLine;
    while (std::getline(file, testCaseLine))
    {
        std::vector<IOPair> tc = parseTestCase(testCaseLine);
        for (size_t i = 0; i < k; ++i) {
            applyTestCase(testCaseLine,tc,testSuiteTree);
        }
    }
}

static void checkSuiteTreeFullyVisited(std::shared_ptr<TestNode> node) {
    if (!node->visited) {
        std::cerr << "Test case not observed in SUT: ";
        printPathToNode(node);
        std::cerr << std::endl;
        std::cerr << "SUT is does not pass the test suite." << std::endl;
        exit(1);
    }

    for (auto& entry : node->nexts) {
        checkSuiteTreeFullyVisited(entry.second);
    }
}

int main(int argc, char* argv[])
{
    parseParameters(argc,argv);

    // create a test suite tree representing all traces in the test suite
    std::shared_ptr<TestNode> testSuiteTree = generateTestSuiteTree();
    testSuiteTree->visited = true;
    
    printTestTree(testSuiteTree);

    // initialize the sut and apply each test case in the test suite
    sut_init();
    applyTestSuite(testSuiteTree);

    printTestTree(testSuiteTree);

    // if the SUT was to be tested w.r.t. language equivalence, it is furthermore necessary
    // to ensure that all behaviours described in the test suite have been observed during
    // testing
    if (conformanceRelation == EQUIVALENCE) {
        checkSuiteTreeFullyVisited(testSuiteTree);
    }
    
    exit(0);
    
}

