#ifndef ECGRAPH_AUTOMATON
#define ECGRAPH_AUTOMATON

#include <vector>
#include <utility>
#include <string>

#include <iostream>

#include <random>

std::default_random_engine getRandomGen(const int custom = -1);


/// @brief Representation of an edge-coloured graph. Can be partiallly filled
class ECGraph{
    /// @brief How many colours are available in the graph
    int nColours;
    /// @brief How many nodes are in the graph
    int nNodes;

    /// @brief All edges saved as adjacency lists. Has size nNodes*nColours
    std::vector<std::vector<int>> edges;

    protected:
        ECGraph();

        void delayedInitalize(int nColours, int nNodes);
        bool isCompatible(int numN, int numC);

        void randomiseGraph(const bool fancy = false, const int sNode=0, 
                            const std::vector<int> acceptingStates = std::vector<int>());

        bool shortRead(std::istream& redF);
        void print(std::ostream& stream, int mode);


    public:
        ECGraph(int nColours, int nNodes);
        ECGraph(std::istream& redF);

        int getNColours();
        int getNNodes();

        int addNode();
        int addNode(std::vector<int> newEdges);

        int adjustEdge(int from, int colour, int to);

        int followEdge(int from, int colour);
        int parse(std::string word, int initialState);
        int parse(const std::vector<int>& word, int initialState);
        int partParse(const std::vector<int>& word, int initialState);

        std::vector<std::pair<int,int>> getUndefinedTransitions();

        void print(std::ostream& stream);
        int read(std::istream& redF, bool constr);

};

class Automaton: public ECGraph{
    
    /// @brief The starting state from which e.g. parsing happens
    int initialState;
    /// @brief For each state, whether it is appecting (True) or rejecting (False)
    std::vector<bool> acceptingStates;

    bool canReachAccepting(int currentState, std::vector<bool>& visited, std::vector<bool>& canReach);
    void completeFancyRandomisation(const bool forceTrueSink=false);
    
    public:
        Automaton(int nStates, int alphabetSize, const int initialState = 0, 
                    const std::vector<int>& acceptingStateList = std::vector<int>{0}, 
                    const bool randomise=true, const bool forceSink=false);
        Automaton(int alphabetSize, std::vector<int>& signature, int sigLength);
        Automaton(Automaton a, Automaton b);
        Automaton(std::istream& redF);

        void print(std::ostream& stream);

        bool isWellDefined();

        void getSinkStatus(std::vector<bool>& wellConnected);

        std::string generateTrace(int length, int style=0);

        std::pair<bool,int> parse(std::string word);

};

#endif