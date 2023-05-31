#ifndef INTERFIX_GRAPH
#define INTERFIX_GRAPH

#include "automaton.h"
#include "unionfind.h"

#include <vector>
#include <unordered_set>
#include <queue>
#include <list>
#include <utility>

#include <string>
#include <sstream>

#include <chrono>

/// @brief An edge coloured graph, which can be merged using the Red-Blue framework
class InterfixGraph: public ECGraph{

    
    /// @brief Whether the graph has been fully merged
    bool collapsed = false;
    /// @brief Whether the graph records frequency information
    bool hasFreq = false;
    /// @brief How much time the merging took
    std::chrono::duration<double> timeToMerge;
    /// @brief Which heuristic is used to evaluate merges
    int scoreMode = 0;
    /// @brief Tracks how many states where merged whose futures are also the same in one iteration
    int futureConsist;
    /// @brief Tracks how many states where merged in one iteration
    int statesMerged;
    /// @brief Tracks the edsm heuristic for the current iteration
    int edsm;
    /// @brief Tracks how many new edges where introduced in the current iteration
    int newEdges;

    bool freqComp(int i, int j);

    bool pairComp(const std::pair<int,int>& a, const std::pair<int,int>& b);

    virtual void delayedInitalize(int nSymbols, const int nStates=1);

    void enrichWithPrefixTrace(std::vector<int>& trace, int posPrefix, const int currentState, const int offset);

    void freqFromTrace(const std::string& trace, int style=0, bool interfix=true, int freqType=2);

    int mapSymb(int symbol, bool prefix);

    std::string dotNodeType(const int node, bool redOnly=false);

    void recursiveDotWrite(std::ostream& os, std::vector<bool>& visited, int cNode, 
                           std::string name, std::vector<std::string>& nameReg,
                           std::queue<std::pair<int,std::string>>& toExplore,
                           bool right, bool left, bool redOnly=false);

    virtual void addRedState(int state);

    void resetScores();

    virtual int getMaxBlue();

    float scoreMerge();

    bool consistent(int u, int v, int strat=0);

    bool merge(int u, int v, UF& me, int strat=0);

    void finalise();

    float probFromTrace(std::vector<int>& trace, std::vector<std::vector<float>>& memory, 
                        int offFront=0, int offEnd=0, int side=-1);

    int partParse(const std::vector<int>& word, int initialState);

    virtual bool checkTermination();

    protected:

        /// @brief How many symbols the graph has. Internally twice as many colours are used (left and right Cayley graph).
        int nSymbols;
        /// @brief What the root state representing the empty word is
        int rootState = 0;
        /// @brief Wheter the structure is ready for merging
        bool finalised = false;
        /// @brief List of all accepting states
        std::vector<int> acceptedStateList;
        /// @brief List of all rejecting states
        std::vector<int> rejectedStateList;

        /// @brief The Union-Find data structure handling the merging
        UF merger;
        /// @brief A vector recording for each state which red-blue framework colour it is at that moment
        std::vector<int> colour;
        /// @brief A list of all red states according to the red-blue framework
        std::list<int> red;

        void enrichWithSuffixTrace(std::vector<int>& trace, const int pos, const int currentState);

        std::pair<int,bool> traceToSeq(const std::string& trace, std::vector<int>& seq, int style=0);

        int getSColour(int state);

    public:
        InterfixGraph(int nSymbols, const int nStates=1);

        InterfixGraph(InterfixGraph& ig);

        InterfixGraph(const std::string& filepath);

        InterfixGraph();

        void printInfo(std::ostream& stream);

        int followEdge(int from, int colour);

        int parse(const std::vector<int>& word, int initialState);

        int addNode();

        void setScoreMode(int mode=0);

        void collapse(int strat=0, bool printDebug=false);

        void enrichWithFileOfTraces(const std::string& filepath, int style=0, bool recFreq=false, bool interfix=true);

        int getAcceptanceStatus(const int node);

        virtual int getAcceptanceStatus(std::vector<int>& trace, bool vote=false);

        int evaluateAcc(const std::string& trace, int style=0, bool vote=false);

        float testAccFile(const std::string& filepath, std::vector<int>& details, int style=0, bool vote=false);

        void getProbDistForFile(const std::string& filepath, std::vector<float>& probs);

        virtual void enrichWithTrace(const std::string& trace, int style=0, bool interfix=true);

        int to_dot_file(const std::string& filepath, bool right=true, bool left=true, bool redOnly=false);

        int to_txt_file(const std::string& filepath);

};

/*
void load_ground_truth(const std::string& filepath, std::vector<float>& dist);

void retest_model(std::string& modelpath, std::string& testfilepath, std::string& outpath, bool vote=false);

void retest_dataset(std::string& dataModelPath, int nFiles, std::string& dataTestpath,
                    std::string& modelName, std::string& outpath, bool vote=false);

float print_test(InterfixGraph& ifg, std::ofstream& os, std::string& testfilepath, bool vote=false);

void train_test_dataset(std::string& in_path, std::string& out_path, int nFiles, int heuristic, 
                        bool aut, bool sm, bool vote=false);

int split_file_into_train_test(const std::string& input, const std::string& out_train, const std::string& out_test);

float perplexity(std::vector<float>& truth, std::vector<float>& predict);

void generateCompositeDataSet(std::string& folderpath, int nFiles, std::vector<int> filesizes);
*/

#endif