#include <vector>
#include <queue>
#include <deque>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <chrono>
#include <random>

#include <functional>
#include <utility>
#include <algorithm>
#include <numeric>

#include "automaton.h"

using std::chrono::system_clock; using std::default_random_engine; using std::uniform_int_distribution;
using std::istream; using std::ostream; using std::istringstream; using std::ostringstream;
using std::vector; using std::pair; using std::string; using std::cout; using std::cin;


/// @brief Returns a random enging, if no seed is given it uses the current system time as seed
/// @param custom Seed for the random enging
/// @return a default random engine
default_random_engine getRandomGen(const int custom){

    unsigned seed;

    if(custom == -1){
        seed = system_clock::now().time_since_epoch().count();
        seed = seed % 10000;
        //cout << seed << "\n";
    }else{
        seed = custom;
    }
    
    default_random_engine gen(seed);

    return gen;
}


/// @brief Initialise an edge coloured graph with no edges
/// @param nColours number of colours the graph has
/// @param nNodes number of nodes the graph has
ECGraph::ECGraph(int nColours, int nNodes){
    this->nColours = nColours;
    this->nNodes = nNodes;

    edges = vector<vector<int>>(nNodes, vector<int>(nColours, -1));
}

/// @brief Initialise edge coloured graph by reading it from a file. Compatible with files produced by write
/// @param redF Input stream pointing to the start of the graph
ECGraph::ECGraph(istream& redF){
    this->read(redF, true);
}

/// @brief Initialise edge coloured graph with no further information. Best used together with delayedInitialise
ECGraph::ECGraph(){
    nColours = -1;
    nNodes = -1;
}

/// @brief Resets the graph to one with the specified number of nodes and colours as well as no edges
/// @param nColours Number of colours of the new graph
/// @param nNodes Number of nodes of the new graph
void ECGraph::delayedInitalize(int nColours, int nNodes){
    this->nColours = nColours;
    this->nNodes = nNodes;

    edges = vector<vector<int>>(nNodes, vector<int>(nColours, -1));
}

/// @brief Getter for nColours
/// @return Number of colours in the graph
int ECGraph::getNColours(){
    return nColours;
}

/// @brief Getter for nNodes
/// @return Number of nodes in the graph
int ECGraph::getNNodes(){
    return nNodes;
}

/// @brief Checkes if the graph has the specififed number of colours and nodes
/// @param numN Number of nodes the graph is compared to
/// @param numC Number of colours the graph is compared to
/// @return True if both values are equal, False otherwise
bool ECGraph::isCompatible(int numN, int numC){
    return numC == nColours && numN == nNodes;
}

/// @brief adds a new node with default edges (going to node -1) for every colour
/// @return the id of the new node
int ECGraph::addNode(){
    nNodes++;
    edges.push_back(vector<int>(nColours, -1));
    return nNodes - 1;
}


/// @brief Add a new node with a predefined set of edges
/// @param newEdges Set of new edges. Length of vector should be equal to nColours
/// @return -1 on failure and the id of the new edge otherwise
int ECGraph::addNode(vector<int> newEdges){
    if((int) newEdges.size() > nColours){
        cout << "Error: Adding new node failed. New node has to many outgoing edges. \n";
        return -1;
    }else if((int) newEdges.size() < nColours){
        cout << "Warning: New node has to few edges. Missing edges will be filled with default edge -1. \n";
        while((int) newEdges.size() < nColours){
            newEdges.push_back(-1);
        }
    }
    edges.push_back(newEdges);
    nNodes++;
    return nNodes - 1;
}

/// @brief Adjusts an edge to have the properties in the argument
/// @param from Node from which the adjusted edge originates
/// @param colour colour of the new adjusted edge
/// @param to Node to which the new adjusted edge points
/// @return 1 on failure and 0 on success
int ECGraph::adjustEdge(int from, int colour, int to){
    if(colour >= nColours){
        cout << "Edge (" << from << ", " << to << ", " 
            << colour << ") rejected, as " << colour << " is not a valid colour \n";
        return 1;
    }
    if(from >= nNodes){
        cout << "Edge (" << from << ", " << to << ", " 
            << colour << ") rejected, as " << from << " is not a valid node \n";
        return 1;
    }
    if(to >= nNodes){
        cout << "Edge (" << from << ", " << to << ", " 
            << colour << ") rejected, as " << to << " is not a valid node \n";
        return 1;
    }
    edges[from][colour] = to;
    return 0;
}

/// @brief Return the node an edge from "from" of colour "colour" points to
/// @param from 
/// @param colour 
/// @return -1 on failure, id of the pointed to node otherwise
int ECGraph::followEdge(int from, int colour){
    if(from < 0 || from >= nNodes){
        cout << "Warning: Undefined node " << from << " was passed to followEdge. \n";
        return -1;
    }else if(colour < 0 || colour >= nColours){
        cout << "Warning: Undefined colour " << colour << " was passed to followEdge. \n";
        return -1;
    }
    return edges[from][colour];
}

/// @brief Identifies and returns all edges which haven't been defined yet
/// @return A vector of all <origin,colour> pairs, which currently don't point to a node
vector<pair<int,int>> ECGraph::getUndefinedTransitions(){
    vector<pair<int,int>> undefinedTransitions = vector<pair<int,int>>();
    for(int i = 0; i < nNodes; i++){
        for(int j = 0; j < nColours; j++){
            if(edges[i][j] == -1){
                undefinedTransitions.push_back(pair<int,int>(i,j));
            }
        }
    }
    return undefinedTransitions;
}

/// @brief Takes a string of colours and an initial state and parses it as if the edge coloured graph was a DFA
/// @param word The sequence of colours to be parsed given as string
/// @param initialState The initial state of the DFA
/// @return -1 if it encounters an undefined transition, the id of the state it is in after parsing the last character otherwise
int ECGraph::parse(std::string word, int initialState){
    int currentState = initialState;
    int currentChar;
    istringstream ss(word);

    do{
        ss >> currentChar;
        currentState = followEdge(currentState, currentChar);
        if(currentState == -1){
            cout << "Warning: Trace led to an undefined state. \n";
            return -1;
        }
    }while(ss);

    return currentState;
}

/// @brief Takes a sequence of colours and an initial state and parses it as if the edge coloured graph was a DFA
/// @param word Sequence of colours given as a vector
/// @param initialState The initial state of the DFA
/// @return -1 if it encounters an undefined transition, the id of the state it is in after parsing the last character otherwise
int ECGraph::parse(const std::vector<int>& word, int initialState){
    int currentState = initialState;
    for(auto it = word.begin(); it != word.end(); ++it){
        currentState = followEdge(currentState, *it);
        if(currentState == -1){
            //cout << "Warning: Trace led to an undefined state. \n";
            return -1;
        }
    }
    return currentState;
}

/// @brief Parses a trace until it encounters the end or an unkown transition.
/// @param word input word/trace that is parsed
/// @param initialState to root state from where the parsing starts
/// @return The last well defined state encountered while parsing
int ECGraph::partParse(const std::vector<int>& word, int initialState){
    int currentState = initialState;
    int nextState;
    for(auto it = word.begin(); it != word.end(); ++it){
        nextState = followEdge(currentState, *it);
        if(nextState == -1){
            return currentState;
        }else{
            currentState = nextState;
        }
    }
    return currentState;
}

/// @brief Create a random edge coloured graph. All edges will be defined if fancy is False, Fancy=true might leave eges undefined
/// @param fancy False - Purely random transitions. True - All nodes a reachable from a starting state
/// @param sNode Node the will act as starting state for fancy randomisation
/// @param acceptingStates Currently not used due to a bug
void ECGraph::randomiseGraph(const bool fancy, const int sNode, 
                             const vector<int> acceptingStates){

    auto gen = getRandomGen();

    if(!fancy){
        uniform_int_distribution<int> dist(0,(nNodes-1));
        auto rEdge = bind(dist, gen);

        for(auto it=edges.begin(); it != edges.end(); ++it){
            for(auto it2=(*it).begin(); it2 != (*it).end(); ++it2){
                *it2 = rEdge();
            }
        }

        return;
    }else{
        uniform_int_distribution<int> dist(1,(nColours));
        auto rnOut = bind(dist, gen);

        vector<int> colours (nColours);
        std::iota(colours.begin(), colours.end(), 0);

        std::deque<int> nodes(nNodes);
        std::queue<int> q;
        
        std::iota(nodes.begin(), nodes.end(), 0);
        nodes.erase(nodes.begin()+sNode);
        std::shuffle(nodes.begin(), nodes.end(), gen);

        q.push(sNode);

        while(!q.empty() && !nodes.empty()){

            std::deque cColours(colours.begin(), colours.end());
            std::shuffle(cColours.begin(), cColours.end(), gen);

            int nOut = rnOut();

            do{

                edges[q.front()][cColours.front()] = nodes.front();
                q.push(nodes.front());
                nodes.pop_front();
                cColours.pop_front();

                nOut--;

            }while(!cColours.empty() && nOut > 0 && !nodes.empty());

            q.pop();

        }
    }
    
}

/// @brief Prints partial or full information about the ECGraph to an output stream
/// @param stream Outstream that will be printed to
/// @param mode 0 -- Prints all information, 1-- only prints edges, 2 -- only prints nNodes and nColours
void ECGraph::print(ostream& stream, int mode){

    if(mode != 1){
        stream << nNodes << " " << nColours << "\n";
    }
    
    if(mode != 2){
        for(int from = 0; from < nNodes; from++){
            for(int c = 0; c < nColours; c++){
                stream << from << " " << c << " " << edges[from][c] << "\n";
            }
        }
    }

}

/// @brief Prints the entire ECGraph to an output stream
/// @param stream 
void ECGraph::print(ostream& stream){
    print(stream, 0);
}

/// @brief Reads graph from an input file. Only reads edges.
/// @param redF Input stream to read from
/// @return True if a full graph was read, False otherwise
bool ECGraph::shortRead(istream& redF){
    string line;

    bool complete =true;

    int from;
    int colour;
    int to;
    for(int nLines = 0; nLines < nNodes * nColours; nLines++){
        if(getline(redF, line)){
            istringstream ss(line);
            ss >> from >> colour >> to;
            edges[from][colour] = to;
        }else{
            complete = false;
            break;
        }
    }

    return complete;
}

/// @brief Makes the ECGraph into a copy of one defined in an input stream
/// @param redF Input stream that is read from
/// @param constr If true, will reject the input, if it doesn't align with the already defined number of nodes and colours
/// @return 0 -- on success, 1 -- on failure
int ECGraph::read(istream& redF, bool constr){
    int numC;
    int numN;
    string line;

    if(getline(redF, line)){
        istringstream ss(line);
        ss >> numN >> numC;
    }else{
        return 1;
    }
    if(!constr){
        if(!isCompatible(numN, numC)){
            cout << "Error: red number of colours or nodes doesn't match with graph. \n";
            return 1;
        }
    }else{
        nColours = numC;
        nNodes = numN;

        edges = vector<vector<int>>(nNodes, vector<int>(nColours, -1));
    }

    shortRead(redF);

    return 0;
}


/// @brief Checks for current and all successor states, if they can reach an accepting state and saves the answer in canReach
/// @param currentState Current state that is analysed
/// @param visited List of all states. Length needs to be equal to nStates and should be initialised with all false
/// @param canReach List of all states. Length needs to be equal to nStates, will contain answer for all successors. 
/// @return If the current state can reach an accepting state
bool Automaton::canReachAccepting(int currentState, vector<bool>& visited, vector<bool>& canReach){

    if(visited[currentState]){
        return canReach[currentState];
    }

    visited[currentState] = true;

    if(canReach[currentState]){
        return true;
    }

    bool can = false;
    int transition = 0;
    while(!can && transition < getNColours()){

        int nextState = followEdge(currentState, transition);

        if(nextState != -1){
            can = can || canReachAccepting(nextState, visited, canReach);
        }
        transition++;
    }
    canReach[currentState] = can;
    return can;
}


/// @brief Initialise Automaton. Either an automaton where all edges are undefined or one, where they are randomised
/// @param nStates Number of states in the automaton
/// @param alphabetSize Size of the alphabet of the automaton
/// @param initialState Root state of the automaton
/// @param acceptingStateList List of all state ids, that will be accepting
/// @param randomise Wether edges should be randomised. If true, all states are reachable from the starting state
/// @param forceSink If true the randomisation will guarantee at least one sink state exists
Automaton::Automaton(int nStates, int alphabetSize, const int initialState, 
            const vector<int>& acceptingStateList, 
            const bool randomise, const bool forceSink) 
            : ECGraph(alphabetSize, nStates){

    this->acceptingStates = vector<bool>(nStates, false);
    for(auto it = acceptingStateList.begin(); it != acceptingStateList.end(); ++it){
        acceptingStates[*it] = true;
    }

    if(initialState < 0 || initialState >= nStates){
        cout << "Error: Invalid initial state " << initialState << " was proposed for model with " 
                << nStates << "total states. \n";
    }
    this->initialState = initialState;

    if(randomise){
        randomiseGraph(true, initialState, acceptingStateList);
        completeFancyRandomisation(forceSink);
    }
}

/// @brief Initialises the automaton from an input stream. Compatible with print()
/// @param redF Input stream to be read from
Automaton::Automaton(istream& redF){
    int numC=0;
    int numN=0;
    int initialState=-1;
    vector<int> acceptingStateList = vector<int>();

    string line;

    if(getline(redF, line)){
        istringstream ss(line);
        ss >> numN >> numC;
    }
    if(getline(redF, line)){
        istringstream ss(line);
        ss >> initialState;
    }
    if(getline(redF, line)){
        int accState;
        istringstream ss(line);
        do{
            ss >> accState;
            acceptingStateList.push_back(accState);
        }while(ss);
    }

    delayedInitalize(numC, numN);

    this->acceptingStates = vector<bool>(numN, false);
    for(auto it = acceptingStateList.begin(); it != acceptingStateList.end(); ++it){
        acceptingStates[*it] = true;
    }

    if(initialState < 0 || initialState >= numN){
        cout << "Error: Invalid initial state " << initialState << " was proposed for model with " 
                << numN << "total states. \n";
        initialState = -1;
    }
    this->initialState = initialState;

    bool complete = shortRead(redF);

    if(!complete){
        cout << "Error: Incomplete state transition list was provided. \n";
    }

}

/// @brief Creates an Automaton which accepts all words of the form "x * y", where a accepts x and b accepts y 
/// @param a First Automaton to be composed
/// @param b Second Automaton to be composed
Automaton::Automaton(Automaton a, Automaton b):ECGraph(){

    int aC = a.getNColours();
    int bC = b.getNColours();
    int aN = a.getNNodes();
    int bN = b.getNNodes();

    int nC = aC + bC;
    int nN = aN + bN + 1;

    delayedInitalize(nC, nN);

    initialState = a.initialState;
    for(int i=0; i < aN; i++){
        for(int c=0; c < aC; c++){
            int next = a.followEdge(i, c);
            adjustEdge(i, c, next);
        }
    }
    for(int i=0; i < bN; i++){
        for(int c=0; c < bC; c++){
            int next = b.followEdge(i, c);
            adjustEdge(aN+i, aC+c, aN+next);
        }
    }
    for(int i=0; i < aN; i++){
        if(a.acceptingStates[i]){
            for(int c=0; c < bC; c++){
                int next = b.followEdge(0, c);
                adjustEdge(i, aC+c, aN+next);
            }
        }
    }
    acceptingStates = std::vector<bool>(nN, false);
    for(int i=0; i < aN; i++){
        acceptingStates[i] = a.acceptingStates[i];
    }
    for(int i=0; i < bN; i++){
        acceptingStates[aN+i] = b.acceptingStates[i];
    }

    for(int i=0; i < nN; i++){
        for(int c=0; c < nC; c++){
            if(followEdge(i, c) == -1){
                adjustEdge(i, c, nN-1);
            }
        }
    }

}

/// @brief Initialises an Automaton that accepts all words of the form "Sigma^* * signature * Sigma^*"
/// @param alphabetSize Size of the alphabet
/// @param signature the sequence that serves as a signature
/// @param sigLength the length of the signature
Automaton::Automaton(int alphabetSize, std::vector<int>& signature, int sigLength):ECGraph(alphabetSize, sigLength+1){

    if(signature.empty()){
        auto gen = getRandomGen();
        uniform_int_distribution<int> dist(0,(alphabetSize-1));
        auto rSymb = bind(dist, gen);
        for(int i=0; i < sigLength; i++){
            signature.push_back(rSymb());
        }
    }

    initialState = 0;
    int inSymb = signature[0];
    int priorLength = 0;
    for(int i=0; i < sigLength; i++){
        if(signature[i] == inSymb){
            priorLength++;
        }
    }

    for(int i = 0; i < sigLength; i++){
        for(int c = 0; c < alphabetSize; c++){
            if(signature[i] == c){
                adjustEdge(i, c, i+1);
            }else{ 
                if(i == priorLength && c == inSymb){
                    adjustEdge(i, c, i);
                }else if(signature[0] == c){
                    adjustEdge(i, c, 1);
                }else{
                    adjustEdge(i, c, 0);
                }
            }
        }
    }
    for(int c = 0; c < alphabetSize; c++){
        adjustEdge(sigLength, c, sigLength);
    }

    acceptingStates = std::vector<bool>(sigLength+1, false);
    acceptingStates[sigLength] = true;

}

/// @brief Print automaton to an output stream
/// @param stream Stream that is printed to
void Automaton::print(ostream& stream){
    ECGraph::print(stream, 2);
    stream << initialState << "\n";
    for(int i =0; i < (int) acceptingStates.size(); i++){
        if(acceptingStates[i]){
            stream << i << " ";
        }
    }
    stream << "\n";
    ECGraph::print(stream, 1);
}

/// @brief Checks wether all transitions are defined
/// @return Returns True if well defined, false otherwise
bool Automaton::isWellDefined(){
    return initialState != -1 && (getUndefinedTransitions().size() == 0);
}


/// @brief Checks for all states, wether they are sinks
/// @param wellConnected should be of length nStates. Will be false for all sinks, true otherwise
void Automaton::getSinkStatus(vector<bool>& wellConnected){
    vector<bool> canReach(acceptingStates);
    vector<bool> visited(getNNodes(), false);

    for(int i = 0; i < getNNodes(); i++){
        canReachAccepting(i, visited, canReach);
    }

    wellConnected = vector<bool>(canReach);

}

/// @brief Helper function for fancy randomisation
/// @param forceTrueSink If true, will contain at least one sink
void Automaton::completeFancyRandomisation(const bool forceTrueSink){

    auto gen = getRandomGen();
    uniform_int_distribution<int> dist(0,(getNNodes()-1));
    auto rEdge = bind(dist, gen);

    int trueSinkCandidate = rEdge();

    if(forceTrueSink){
        bool foundSink;
        for(int i = 0; i < getNNodes(); i++){
            bool flag = false;
            for(int j = 0; j < getNColours(); j++){
                flag = flag || followEdge(i, j) != -1;
            }
            if(!flag){
                for(int j = 0; j < getNColours(); j++){
                    adjustEdge(i, j, i);
                }
                foundSink = true;
                trueSinkCandidate = i;
                break;
            }
        }
        if(!foundSink){
            cout << "Warning: Didn't find valid candidate for true sink. \n";
        }
    }

    vector<bool> wellConnected(getNNodes());

    //connect sinks to accepting states
    std::deque<int> sinks;
    std::deque<int> sources;
    do{

        getSinkStatus(wellConnected);

        sinks.clear();
        sources.clear();

        for(int i=0; i < getNNodes(); i++){
            if(wellConnected[i]){
                sources.push_back(i);
            }else{
                if(i != trueSinkCandidate){
                    sinks.push_back(i);
                }
            }
        }

        std::shuffle(sources.begin(), sources.end(), gen);
        std::shuffle(sinks.begin(), sinks.end(), gen);

        bool flag = false;
        while(!sinks.empty() && !flag){
            for(int i=0; i < getNColours(); i++){
                if(followEdge(sinks.front(), i) == -1){
                    adjustEdge(sinks.front(), i, sources.front());
                    flag = true;
                    break;
                }
            }
            sinks.pop_front();
        }

    }while(!sinks.empty() && (int) sources.size() < getNNodes()-1);

    auto toDetermine = getUndefinedTransitions();
    //randomly assign the remaining edges
    while(!toDetermine.empty()){
        adjustEdge(toDetermine.back().first, toDetermine.back().second, rEdge());
        toDetermine.pop_back();
    }

}

/// @brief Returns a trace that can be generated by the automaton
/// @param length Requested length of the generated trace
/// @param style 0 -- length, accepted/rejected, trace; 1 -- accepted/rejected, length, trace
/// @return 
string Automaton::generateTrace(int length, int style){
    ostringstream oss;
    int currentChar;

    auto gen = getRandomGen();
    uniform_int_distribution<int> dist(0,(getNColours()-1));
    auto rChar = bind(dist, gen);

    int currentState = initialState;
    oss << length << " ";

    for(int i=length; i > 0; i--){
        currentChar = rChar();
        oss << currentChar << " ";
        currentState = followEdge(currentState, currentChar);
    }

    oss << "\n";

    string trace = oss.str();
    string accS;

    switch(style){

        case 0:
            accS = acceptingStates[currentState]?" 1" :" 0";
            trace.insert(std::to_string(length).size(), accS);
            break;

        case 1:
            accS = acceptingStates[currentState]?"1 " :"0 ";
            trace.insert(0, accS);
            break;

        default:
            accS = acceptingStates[currentState]?" 1" :" 0";
            trace.insert(std::to_string(length).size(), accS);
            break;

    }
        

    return trace;
}

/// @brief Tries to parse a string as Automaton
/// @param word String that is to be parsed
/// @return A pair of wether the word is accepted by the automaton and the state in which the trace ended. Will be False,-1 if it encounters an undefined transition
pair<bool,int> Automaton::parse(string word){
    int finalState = ECGraph::parse(word, initialState);

    if(finalState == -1){
        return pair<bool,int>(false, -1);
    }

    bool success = acceptingStates[finalState];

    return pair<bool,int>(success,finalState);
}