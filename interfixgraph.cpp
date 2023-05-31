#include "automaton.h"
#include "interfixgraph.h"
#include "unionfind.h"

#include <vector>
#include <string>
#include <unordered_set>
#include <set>
#include <queue>
#include <iterator>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <algorithm>
#include <functional>
#include <utility>
#include <chrono>

#include <cmath>
#include <numeric>
#include <random>


/// @brief Maps alphabet ids to the internally used one depending on if they're used as pre- or appending symbols
/// @param symbol the alphabet id to be mapped
/// @param prefix If true, symbol is in the left Cayley graph, if false in the right Cayley graph
/// @return the mapped id
int InterfixGraph::mapSymb(int symbol, bool prefix){
    if(symbol >= nSymbols){
        return prefix? symbol: symbol % nSymbols;
    }else{
        return prefix? symbol + nSymbols: symbol;
    }
}

/// @brief Initialise an Interfixgraph with a set alphabet size and number of states but no edges
/// @param nSymbols Size of the alphabet
/// @param nStates Number of states
InterfixGraph::InterfixGraph(int nSymbols, const int nStates): ECGraph(nSymbols * 2, nStates){
    rootState = 0;
    this->nSymbols = nSymbols;
    merger = UF(nStates, nSymbols * 2);
}

/// @brief Copy constructor
/// @param ig Interfix Graph that is copied
InterfixGraph::InterfixGraph(InterfixGraph& ig):ECGraph(){

    ECGraph::delayedInitalize(ig.getNColours(), ig.getNNodes());
    merger = UF(getNNodes(), getNColours());
    nSymbols = ig.getNColours() / 2;
    rootState = ig.rootState;
    for(int i = 0; i < getNColours(); i++){
        int acc = ig.getAcceptanceStatus(i);
        if(acc != -1){
            i ? acceptedStateList.push_back(i) : rejectedStateList.push_back(i);
        }
        merger.setAcc(i, acc);
        for(int c = 0; c < getNColours(); c++){
            adjustEdge(i, c, ig.followEdge(i,c));
        }
    }
    finalised = true;

}

/// @brief File read constructor. Compatible with to_txt_file()
/// @param filepath path to the file from which the Interfix Graph should be read
InterfixGraph::InterfixGraph(const std::string& filepath):ECGraph(){

    std::ifstream fs(filepath);

    if(fs.is_open()){

        int nNodes;
        int nSymbs;

        std::string line;

        if(std::getline(fs, line)){
            std::istringstream ss(line);
            ss >> nNodes >> nSymbs;
            ECGraph::delayedInitalize(2*nSymbs, nNodes);
            this->nSymbols = nSymbs;
            merger = UF(getNNodes(), getNColours());
            merger.initaliseAcc(getNNodes());
            merger.initaliseFreq(getNNodes(), getNColours());
            this->rootState = 0;
            for(int i=0; i < getNNodes(); i++){
                if(std::getline(fs, line)){
                    int current;
                    int acc;
                    float freq;
                    ss = std::istringstream(line);
                    ss >> current >> acc >> freq;
                    merger.setAcc(current, acc);
                    merger.incFreq(current, freq);
                    for(int c=0; c < getNColours(); c++){
                        if(std::getline(fs, line)){
                            ss = std::istringstream(line);
                            int colour;
                            int next;
                            ss >> current >> colour >> next >> freq;
                            adjustEdge(current, colour, next);
                            merger.incOutFreq(current, colour, freq);
                        }else{
                            std::cout << "Error: couldn't read information about edge " << i << ", " << c
                                      << " from file " << filepath << "\n";
                            return;
                        }
                    }
                }else{
                    std::cout << "Error: couldn't read information about node " << i << " from file " << filepath << "\n";
                    return;
                }
            }
        }else{
            std::cout << "Error: couldn't read first line from input file " << filepath << "\n";
            return;
        }
    }else{
        std::cout << "Error: couldn't read from input file " << filepath << "\n";
        return;
    }
    finalised = true;
}

/// @brief Generic constructor. Needs call of delayed initialise before being useable
InterfixGraph::InterfixGraph():ECGraph(){
    rootState = -1;
}

/// @brief Print the summarised info of state merging. Will not work properly if collapse() hasn't been called yet
/// @param stream Output stream that is printed to
void InterfixGraph::printInfo(std::ostream& stream){
    if(!collapsed){
        stream << "A problem occured during learning. \n";
        return;
    }
    stream << std::fixed << std::setprecision(2);
    stream << "\nInitial number of states: " << getNNodes() << "\n";
    stream << "Number of states after merging: " << merger.getNSets() << "\n";
    stream << "Time to collapse: " << timeToMerge.count() << "s\n";
}

/// @brief Updates structures so the interfix graph works with the given alphabet size and number of states
/// @param nSymbols Alphabet size
/// @param nStates Number of states
void InterfixGraph::delayedInitalize(int nSymbols, const int nStates){
    ECGraph::delayedInitalize(2*nSymbols, nStates);
    rootState = 0;
    this->nSymbols = nSymbols;
    merger = UF(nStates, nSymbols*2);
}

/// @brief Less then comparator over the frequency of states
/// @param i left-hand side of comparison
/// @param j right-hand side of comparison
/// @return True if frequency of i is smaller then frequency of j. False otherwise
bool InterfixGraph::freqComp(int i, int j){
    return merger.getFreq(i) < merger.getFreq(j);
}

/// @brief Updates structures, such that the Interfix Graph has an additional state
/// @return Id of the new state
int InterfixGraph::addNode(){
    merger.addSet();
    return ECGraph::addNode();
}

/// @brief Returns the state the given edge points to. Compatible with merging.
/// @param from Node from which the edge originates
/// @param colour Colour/symbol of the edge
/// @return -1 if edge doesn't exist and id of the state the edge points to otherwise
int InterfixGraph::followEdge(int from, int colour){
    int next = ECGraph::followEdge(from, colour);
    if(next == -1){
        int mState = merger.getMergedEdge(from, colour);
        if(mState != -1){
            return followEdge(mState, colour);
        }else{
            return -1;
        }
    }else{
        return merger.find(next);
    }
}

/// @brief Parses a word through the right graph
/// @param word Word given as a sequence of symbols
/// @param initialState State from which the parsing starts
/// @return -1 if an udefined transition is encountered, otherwise the id of the state in which the parsing terminates
int InterfixGraph::parse(const std::vector<int>& word, int initialState){
    int currentState = initialState;
    for(auto it = word.begin(); it != word.end(); ++it){
        currentState = followEdge(currentState, *it);
        if(currentState == -1){
            return -1;
        }
    }
    return currentState;
}

/// @brief Parses a word until it encounters an undefined transition
/// @param word Word given as a sequence of symbols
/// @param initialState state from which the parsing starts
/// @return The id of the last defined state encountered while parsing the word
int InterfixGraph::partParse(const std::vector<int>& word, int initialState){
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

/// @brief Returns the colour (Red, Blue, White) of a state according to the Red-Blue merging framework
/// @param state Id of the state whose colour is queried
/// @return -1 -- White, 0 -- Blue, 1 -- Red
int InterfixGraph::getSColour(int state){
    if(colour.empty()){
        colour = std::vector<int>(getNNodes(), -1);
    }
    return colour[merger.find(state)];
}

/// @brief Identify the blue state with the highest frequency
/// @return The id of the blue state with the highest frequency
int InterfixGraph::getMaxBlue(){
    bool hasFreq = merger.hasFreq();
    int i = -1;
    int max = -1;
    for(auto it = red.begin(); it != red.end(); ++it){
        i = merger.find(*it);
        for(int c = 0; c < getNColours(); c++){
            int next = followEdge(i, c);
            if(next != -1 && getSColour(next) != 1){
                if(!hasFreq){
                    return next;
                }else{
                    if(max == -1){
                        max = next;
                    }else{
                        if(merger.getFreq(next) > merger.getFreq(max)){
                            max = next;
                        }
                    }
                }
            }
        }
    }

    return max;
}

/// @brief Change the colour of a given state to red. Updates list of red states and the colour of all successors.
/// @param state Id of the state that will be turned into a red state
void InterfixGraph::addRedState(int state){
    if(colour.empty()){
        colour = std::vector<int>(getNNodes(), -1);
    }
    int i = merger.find(state);
    colour[i] = 1;
    //Insert state into correct ordered position
    bool biggest = true;
    for(auto it = red.begin(); it != red.end(); ++it){
        if(merger.getFreq(i) <= merger.getFreq(*it)){
            red.insert(it, i);
            biggest = false;
            break;
        }
    }
    if(biggest){
        red.insert(red.end(), i);
    }
    for(int c = 0; c < getNColours(); c++){
        int next = followEdge(i, c);
        if(next != -1){
            if(colour[next] == -1){
                colour[next] = 0;
            }
        }
    }
}

/// @brief Resets the scores used for metrics that evaluate the quality of a merge
void InterfixGraph::resetScores(){

    futureConsist = 0;
    statesMerged = 0;
    edsm = 0;
    newEdges = 0;

}

/// @brief Checks whether all states are red
/// @return true if all states are red, false otherwise
bool InterfixGraph::checkTermination(){
    return merger.getNSets() > (int) red.size();
}


/// @brief Performs a Red-Blue-Framework state-merging algorithm on the structure. Works in place.
/// @param strat 0 -- default with positive and negative traces; 1 -- Alergia94 with only positive traces; 2 -- doesn't merge states who have all positive or negative futures
/// @param printDebug If true, will generate debug information in debug folder and stop after 150 iterations
void InterfixGraph::collapse(int strat, bool printDebug){

    //Start measuring the time
    auto start = std::chrono::system_clock::now();

    //Add root as initial red state
    addRedState(rootState);
    //Finalise acceptance status before starting algorithm
    finalise();

    std::string debugfilepath = ".\\debug\\most_recent.txt";
    std::string debugDotpath = ".\\debug\\most_recent_dot_files\\";
    std::ofstream dos(debugfilepath);

    int its = 0;
    int maxits = getNNodes();

    while(merger.getNSets() > (int) red.size() && its <= maxits){

        //Pick blue state with highest frequency
        int b = getMaxBlue();
        if(b == -1){
            std::cout << "Warning: No valid blue states left. \n";
            break;
        }

        float maxScore = -1.0;
        int posMerges = 0;
        std::pair<int,int> bestMerge;

        bool mergeable = false;
        for(auto itr = red.rbegin(); itr != red.rend(); ++itr){
            resetScores();
            bool success = false;
            //If all futures disagree don't merge
            if(strat == 2){
                if(((bool)merger.getNFutureAcc(*itr) != (bool)merger.getNFutureAcc(b)) ||
                ((bool)merger.getNFutureRej(*itr) != (bool)merger.getNFutureRej(b))){
                    success = false;
                }else{
                    success = merge(*itr, b, merger, strat);
                }
            }else{
                success = merge(*itr, b, merger, strat);
            }
            if(success){
                mergeable = true;
                posMerges++;
                if(scoreMode == 0){
                    break;
                }else{
                    float score = scoreMerge();
                    if(score > maxScore){
                        maxScore = score;
                        bestMerge = std::pair<int,int>(*itr, b);
                    }
                }
                
            }
            //Undo merge, as it was either unsuccessful or there are more merges to explore
            merger.undoMergeHist();

        }
        if(mergeable){
            float fr = merger.getFreq(bestMerge.first);
            float fb = merger.getFreq(bestMerge.second);
            resetScores();
            //Accept best merge
            merge(bestMerge.first, bestMerge.second, merger, strat);
            merger.acceptMergeHist();
            if(printDebug){
                dos << "Merge - edsm: " << edsm << "; states merged: " << statesMerged
                    << "; new edges: " << newEdges
                    << "; possible merges " << posMerges
                    << "; fr: " << fr 
                    << "; fb: " << fb << ";\n";
            }
            //resort red states, as frequencies might have changed
            red.sort([this] (int i, int j){ return freqComp(i, j); });
        }else{
            //Add the blue state to red states, as it couldn't be merged with any red state
            addRedState(b);
            if(printDebug){
                dos << "New red - freq: " << merger.getFreq(b) << "\n";
            }
        }

        its++;
        if(printDebug){
            std::string filename = std::to_string(its) + "_its_identified.dot";
            std::string outpath = debugDotpath + filename;
            to_dot_file(outpath, true, true, true);
            if(its == 150){
                return;
            }
        }
        if(its % 10 == 0){
            std::cout << "|";
        }
        if(its % 1000 == 0){
            std::cout << "\n";
            std::cout << "Progress report: \n";
            std::cout << "Iteration: " << its << "\n";
            std::cout << "Number of states: " << merger.getNSets() << "\n";
            std::cout << "Number of red states: " << red.size() << "\n\n";
        }
    }

    if(checkTermination()){
        std::cout << "Erronous state merging occured. Results are faulty. \n"; 
        return;
    }

    //Finish time measurements
    collapsed = true;
    auto end = std::chrono::system_clock::now();
    timeToMerge = end - start;
    dos.close();
}

/// @brief Checks wether two states u,v are consistent for merging.
/// @param u 
/// @param v 
/// @param strat default -- Does accepting status match?; 1 -- Alergia94 check 
/// @return True if u and v are considered consistent for merging, False otherwise
bool InterfixGraph::consistent(int u, int v, int strat){
    int accu = merger.getAcc(u);
    int accv = merger.getAcc(v);
    //Check if merge is consistent
    if(accu + accv == 1){
        return false;
    }
    //Each pair of red states was already determined to be an inconsistent merge.
    //Accordingly the procedure can be stopped early if such a pair is encountered.
    if(getSColour(u) == 1 && getSColour(v) == 1){
        return false;
    }

    //Alergia94 consistency check
    if(strat == 1){
        double alpha = 0.001;
        double frequ = merger.getFreq(u);
        double freqv = merger.getFreq(v);
        double bound = std::sqrt(0.5 * std::log(2.0/alpha)) * (1/std::sqrt(frequ) + 1/std::sqrt(freqv));
        for(int c = 0; c < getNColours(); c++){
            double gamma = (((double) merger.getOutFreq(u,c))/frequ) -  (((double) merger.getOutFreq(v,c))/freqv);
            if(std::abs(gamma) > bound){
                return false;
            }
        }
    }

    return true;
}

/// @brief Recursive merging that checks for the current merge and all merges encountered during determinisation if they're consistent
/// @param u State that is merged
/// @param v Other state that is merger
/// @param me the merger UnionFind structure
/// @param strat Strategy for checking consistency. 0 -- accept/reject; 1 -- Alergia94
/// @return True, if the current and all merges caused due to determinisation are consistent. False otherwise.
bool InterfixGraph::merge(int u, int v, UF& me, int strat){
    int ru = me.find(u);
    int rv = me.find(v);
    int accu = me.getAcc(ru);
    int accv = me.getAcc(rv);
    //Check if merge is consistent
    if(!consistent(ru, rv, strat)){
        return false;
    }
    //Go up in recursion if states are already merged
    if(ru == rv){
        return true;
    }
    //Always merge with the following hierarchy: Red <- Blue <- White
    if(getSColour(rv) > getSColour(ru)){
        std::swap(ru,rv);
        std::swap(accu,accv);
    }
    //Merge
    me.merge(ru,rv);
    me.setAcc(ru, std::max(accu,accv));
    //update heuristic scores
    statesMerged++;
    if(!((bool)me.getNFutureAcc(ru) ^ (bool)me.getNFutureAcc(rv))
        && !((bool)me.getNFutureRej(ru) ^ (bool)me.getNFutureRej(rv))){
            futureConsist++;
    }
    if(accu == accv && accu != -1){
        edsm++;
    }
    //Check for determinism, and update edges for merging
    std::vector<std::pair<int,int>> toDet;
    for(int c = 0; c < getNColours(); c++){
        int nexu = followEdge(ru, c);
        int nexv = followEdge(rv, c);
        //Update edges
        if(nexu == -1 && nexv != -1){
            newEdges++;
            me.mergeEdge(ru, c, rv);
        }
        if(nexu != -1 && nexv == -1){
            newEdges++;
        }
        //Check for non-determinism
        if(nexu != nexv && nexu != -1 && nexv != -1){
            toDet.push_back(std::pair<int,int>(nexu,nexv));
        }
    }
    //Recursive call for determinisation
    bool det = true;
    for(auto it=toDet.begin(); it != toDet.end(); ++it){
        det = det && merge((*it).first, (*it).second, me);
        if(!det){
            return false;
        }
    }
    //If no determinisation is necessary, merge was successful
    return det;
}

/// @brief Set how scores are merged. Default is highest frequnecy states are merged first
/// @param mode 0 -- Highest Frequency first
///             1 -- Highest "future" consistency
///             2 -- Most nodes merged in determinisation
///             3 -- EDSM
/// 	        4 -- Least new edges
///             5 -- Most noes merged, but corrected to minimise new edges
void InterfixGraph::setScoreMode(int mode){
    scoreMode = mode;
}

/// @brief Assigns a score to the merge according to score mode. Higher scores are better.
/// @return The score set in score mode
float InterfixGraph::scoreMerge(){

    switch (scoreMode)
    {
    case 0:
        return 1.0;
        break;

    case 1:
        return ((float) futureConsist) / ((float) statesMerged);
        break;

    case 2:
        return ((float) statesMerged);

    case 3:
        return edsm;

    case 4:
        return getNNodes()*getNColours() - newEdges;

    case 5:
        return getNNodes()*getNColours() - newEdges + statesMerged;
    
    default:
    return 1.0;
        break;
    }

}

/// @brief Extends the graph with information from a left to right trace
/// @param trace the trace which is used for extending
/// @param pos variable for recursion. Current position in the trace
/// @param currentState variable for recursion. Current state in the graph. Initialise with root state
void InterfixGraph::enrichWithSuffixTrace(std::vector<int>& trace, const int pos, const int currentState){
    if(pos >= (int) trace.size()){
        return;
    }

    int nextState = followEdge(currentState, mapSymb(trace[pos], false));
    if(nextState == -1){
        nextState = addNode();
        adjustEdge(currentState, mapSymb(trace[pos], false), nextState);
    }
    enrichWithSuffixTrace(trace, pos+1, nextState);

}

/// @brief Extends the graph with information from a right to left trace
/// @param trace the trace which is used for extending
/// @param posPrefix variable for recursion. Current position in the trace
/// @param currentState variable for recursion. Current state in the graph. Initialise with root state
/// @param offset Cutoff point starting from which the trace is considered
void InterfixGraph::enrichWithPrefixTrace(std::vector<int>& trace, int posPrefix, const int currentState,
                                          const int offset){
    if(posPrefix < 0 || trace.begin() + posPrefix == trace.end()){
        return;
    }
    int nextState = ECGraph::parse(std::vector<int>(trace.begin()+posPrefix, trace.end()-offset), rootState);
    ECGraph::adjustEdge(currentState, mapSymb(trace[posPrefix], true), nextState);
    enrichWithPrefixTrace(trace, posPrefix-1, nextState, offset);
}

/// @brief Convert a string into a sequence written into seq.
/// @param trace String that will be converted
/// @param seq Vector container for the resulting sequence
/// @param style How length and membership information is encoded in the string
/// @return A pair of length,membership for the trace
std::pair<int,bool> InterfixGraph::traceToSeq(const std::string& trace, std::vector<int>& seq, int style){
    std::istringstream ss(trace);
    int length;
    bool accepted;

    switch(style){
        case 0:
            ss >> length >> accepted;
            break;
        case 1:
            ss >> accepted >> length;
            break;
        default:
            ss >> length >> accepted;
    }
    

    seq = std::vector<int>(length);

    for(int i=0; i < length; i++){
        ss >> seq[i];
    }

    return std::pair<int,bool>(length, accepted);
}

/// @brief Extends the graph with the information from a trace
/// @param trace The trace used for extending the graph
/// @param style How the trace is encoded. 0 -- length membership trace, 1 -- membership length trace
/// @param interfix Wether the information from the reverse trace is also recorded (in a prepending fashion)
void InterfixGraph::enrichWithTrace(const std::string& trace, int style, bool interfix){
    finalised = false;
    std::vector<int> seq;
    auto lenAcc = traceToSeq(trace, seq, style);
    int length = lenAcc.first;
    bool accepted = lenAcc.second;

    for(int i=0; i < length; i++){
        enrichWithSuffixTrace(seq, i, rootState);
        if(!interfix){
            break;
        }
    }
    if(interfix){
        for(int i=length-1; i >= 0; i--){
            enrichWithPrefixTrace(seq, i, rootState, (length-1)-i);
        }
    }

    int finalState = ECGraph::parse(seq, rootState);
    accepted? acceptedStateList.push_back(finalState) : rejectedStateList.push_back(finalState);
}


/// @brief Comparator object for pairs of int. Puts them in lexicographical order
struct PairCmp{
    bool operator()(const std::pair<int,int>& a, const std::pair<int,int>& b) const{
        if(a.first < b.first){
            return true;
        }else{
            if(a.first == b.first){
                return a.second < b.second;
            }else{
                return false;
            }
        }
    }
};


/// @brief Given a trace which was already inserted into the APTA or IG
/// this function will add frequency information for the trace.
/// @param trace the trace for which frequency is added
/// @param style the encoding style of the trace. Default (1) is Abbadingo
/// @param interfix True if the structure is an IG, false if an APTA
/// @param freqType Should always be 1 and defaults to 1. 0 for backwards compatibility
void InterfixGraph::freqFromTrace(const std::string& trace, int style, bool interfix, int freqType){
    std::vector<int> seq;
    auto lenAcc = traceToSeq(trace, seq, style);
    int length = lenAcc.first;

    switch(freqType){
        //Count each possible path by increasing its frequency by 1. Won't work with probabilities
        case 0:
            {
            std::unordered_set<int> visitedStates;
            std::set<std::pair<int,int>,PairCmp> visitedEdges;

            visitedStates.insert(rootState);

            //get states and edges visited by appending symbols
            for(int i=0; i < length; i++){
                int cState = rootState;
                int nextState = -1;
                for(int pos=i; pos < length; pos++){
                    nextState = followEdge(cState, mapSymb(seq[pos], false));
                    visitedStates.insert(nextState);
                    visitedEdges.insert(std::pair<int,int>(cState, mapSymb(seq[pos], false)));
                    cState = nextState;
                }
                if(!interfix){
                    break;
                }
            }
            if(interfix){
                //get states and edges visited by prepending symbols
                for(int i=length-1; i >= 0; i--){
                    int cState = rootState;
                    int nextState = -1;
                    for(int pos=i; pos >= 0; pos--){
                        nextState = followEdge(cState, mapSymb(seq[pos], true));
                        visitedStates.insert(nextState);
                        visitedEdges.insert(std::pair<int,int>(cState, mapSymb(seq[pos], true)));
                        cState = nextState;
                    }
                }
            }


            //update frequencies
            for(auto it = visitedStates.begin(); it != visitedStates.end(); ++it){
                merger.incFreq(*it);
            }
            for(auto it = visitedEdges.begin(); it != visitedEdges.end(); ++it){
                merger.incOutFreq((*it).first, (*it).second);
            }
            break;
        }
        //Fractional frequencies. Work as probabilities
        case 1:
        default:

            std::vector<std::pair<int,std::vector<int>>> traceGraph;

            traceGraph.push_back(std::pair<int,std::vector<int>>(rootState, std::vector<int>(getNColours(), -1)));

            //get states and edges visited by appending symbols
            for(int i=0; i < length; i++){
                int cState = 0;
                int nextState = -1;
                for(int pos=i; pos < length; pos++){
                    int c = mapSymb(seq[pos], false);
                    if(traceGraph[cState].second[c] != -1){
                        nextState = traceGraph[cState].second[c];
                    }else{
                        nextState = traceGraph.size();
                        traceGraph.push_back(std::pair<int,std::vector<int>>(followEdge(traceGraph[cState].first, c), 
                                                std::vector<int>(getNColours(), -1)));
                    }

                    traceGraph[cState].second[c] = nextState;
                    cState = nextState;

                }
                if(!interfix){
                    break;
                }
            }
            if(interfix){
                //get states and edges visited by prepending symbols
                for(int i=length-1; i >= 0; i--){
                    int cState = 0;
                    int nextState = 0;
                    for(int pos=i; pos >= 0; pos--){
                        for(int fpos = pos; fpos <= i; fpos++){
                            int c = mapSymb(seq[fpos], false);
                            nextState = traceGraph[nextState].second[c];
                        }
                        int c = mapSymb(seq[pos], true);

                        traceGraph[cState].second[c] = nextState;
                        cState = nextState;
                        nextState = 0;
                    }
                }
            }


            //update frequencies
            /* Frequencies will work as expected when parsing into a prefix tree: Each edge and node on a path gets its
            frequency inreeased by one. For interfix graphs we work with partial frequencies. Each possible path a 
            trace could take through an interfix graph is treated as equally likely. As such only the root and the
            final node of such a trace will get their frequency increased by 1, while all others only get fractional
            increases.
            */
            std::queue<int> states;
            std::vector<float> freqs(traceGraph.size(), 0.0);
            std::vector<bool> visited(traceGraph.size(), false);
            states.push(0);
            bool acc = lenAcc.second;
            freqs[0] = 1.0;
            while(!states.empty()){
                int current = states.front();
                if(!visited[current]){
                    float sucCount = 0.0;
                    if(!interfix){
                        sucCount = 1.0;
                    }else{
                        for(int c=0; c < getNColours(); c++){
                            if(traceGraph[current].second[c] != -1){
                                sucCount += 1.0;
                            }
                        }
                    }
                    for(int c=0; c < getNColours(); c++){
                        int next = traceGraph[current].second[c];
                        if(next != -1){
                            states.push(next);
                            float edgeFreq = freqs[current] / sucCount;
                            freqs[next] += edgeFreq;
                            merger.incOutFreq(traceGraph[current].first, c, edgeFreq);
                            if(!interfix){
                                break;
                            }
                        }
                    }
                    merger.incFreq(traceGraph[current].first, freqs[current]);
                    merger.incFutureAcc(traceGraph[current].first, acc);
                    visited[current] = true;
                }
                states.pop();
            }
            break;
    }
    
}

/// @brief Extends the graph with the information present in all traces in the given file
/// @param filepath Path to the file which contains traces
/// @param style Style in which the traces are encoded. 1 -- is scientific standard
/// @param recFreq If set to true the frequency certain prefixes/interfixes appear will also be recorded
/// @param interfix If false will build a prefix tree, if true will build an interfix graph
void InterfixGraph::enrichWithFileOfTraces(const std::string& filepath, int style, bool recFreq, bool interfix){

    std::ifstream fs (filepath);

    //Consume the lines describing the generating automaton, if they exist
    if(style == 0){
        Automaton a(fs);
        delayedInitalize(a.getNColours());
    }

    std::string trace;
    std::vector<std::string> traces;
    int nTraces;

    switch(style){
        case 0:
            fs >> nTraces;
            break;
        case 1:
            int nSymbols;
            fs >> nTraces >> nSymbols;
            delayedInitalize(nSymbols);
            break;
        default:
            fs >> nTraces;
            break;
    }
    std::getline(fs, trace);

    for(int i = 0; i < nTraces; i++){
        std::getline(fs, trace);
        traces.push_back(trace);
        enrichWithTrace(trace, style, interfix);
    }

    fs.close();

    if(recFreq){
        hasFreq = true;
        merger.initaliseFreq(getNNodes(), getNColours());
        for(int i = 0; i < nTraces; i++){
            freqFromTrace(traces[i], style, interfix);
        }
    }

}

/// @brief Predicts how probable it is that a given trace is generated by the model
/// @param trace Trace, given as sequence of symbols, for which the probability is predicted
/// @param memory For recursion needs to be initialised with length equal to trace
/// @param offFront For recursion, leave default argument
/// @param offEnd For recursion, leave default argument
/// @param side For recursion, leave default argument
/// @return Probability of the trace or 1.0e-30, if impossible (since 0 can cause computational problems)
float InterfixGraph::probFromTrace(std::vector<int>& trace, std::vector<std::vector<float>>& memory, 
                                    int offFront, int offEnd, int side){

    //Probability 0 causes issues in further computation, so a minimal probability is set for traces
    float minProb = 1.0e-30;

    int current = parse(std::vector<int>(trace.begin()+offFront, trace.end()-offEnd), rootState);
    //If an undefined transition is encountered, it is treated, as an impossible trace
    if(current == -1){
        return minProb;
    }

    float prob;

    if(side == 0){
        prob = merger.getOutFreq(current, mapSymb(trace[offFront-1], true)) / merger.getFreq(current);
    }else if(side == 1){
        prob = merger.getOutFreq(current, mapSymb(trace[trace.size()-offEnd], false)) / merger.getFreq(current);
    }else{
        //Compute probability, that trace terminates in this state
        float outFreq = 0.0;
        float totFreq = merger.getFreq(current);
        for(int c = 0; c < getNColours(); c++){
            outFreq += merger.getOutFreq(current, c);
        }
        prob = (totFreq - outFreq) / totFreq;
    }

    if(prob == 0.0){
        return minProb;
    }


    if(trace.size() -offFront -offEnd == 0){
        return prob;
    }else{
        if(memory[offFront][offEnd] >= 0.0){
            return memory[offFront][offEnd];
        }
        float total = prob * (probFromTrace(trace, memory, offFront+1, offEnd, 0) 
                                + probFromTrace(trace, memory, offFront, offEnd+1, 1));
        if(memory[offFront][offEnd] < 0.0){
            memory[offFront][offEnd] = total;
        }
        return total;
    }

}

/// @brief Predicts for a file of traces for each trace the probability, that it is produced by the model. Probabilities are normalised so the results sum to 1
/// @param filepath Path to the file with traces
/// @param probs Vector in which the resulting probabilities are stored
void InterfixGraph::getProbDistForFile(const std::string& filepath, std::vector<float>& probs){

    //minimal probability value, to avoid floating point errors
    float minProb = 1.0e-30;

    std::ifstream fs(filepath);

    if(fs.is_open()){

        std::string line;

        if(std::getline(fs, line)){

            int nTraces, nSymbols;
            std::istringstream ss(line);

            ss >> nTraces >> nSymbols;

            probs = std::vector<float>(nTraces, 0.0);

            for(int i=0; i < nTraces; i++){
                if(std::getline(fs, line)){
                    std::vector<int> seq;
                    auto lenacc = traceToSeq(line, seq, 1);
                    std::vector<std::vector<float>> memory(lenacc.first, std::vector<float>(lenacc.first, -1.0));
                    probs[i] = std::max(probFromTrace(seq, memory), minProb);
                }else{
                    std::cout << "Error: couldn't read trace " << i << " from file at " << filepath << "\n";
                }
            }
            std::cout << "\n";

            float sum = std::accumulate(probs.begin(), probs.end(), 0.0f);
            for(auto it=probs.begin(); it != probs.end(); ++it){
                *it = *it / sum;
            }

        }else{
            std::cout << "Error: couldn't read header from file at " << filepath << "\n";
        }

    }else{
        std::cout << "Error: couldn't read from file at " << filepath << "\n";
    }

}

/// @brief Initialise the merger to be ready for merging
void InterfixGraph::finalise(){

    finalised = true;
    merger.finaliseAcc(acceptedStateList, rejectedStateList);

}

/// @brief Returns wether a given node is accpeting or rejecting
/// @param node 
/// @return 1 -- Nodes is accepting; 0 -- Node is rejecting; -1 -- Status of node is unknown
int InterfixGraph::getAcceptanceStatus(const int node){
    if(!finalised){
        finalise();
    }
    return merger.getAcc(node);
}


/// @brief Returns wether a given node is accepting or rejecting.
/// @param trace 
/// @param vote If vote is yes and the status of a node is unknown its status is the majority class of its successors
/// @return 1 -- Nodes is accepting; 0 -- Node is rejecting; -1 -- Status of node is unknown and vote is false
int InterfixGraph::getAcceptanceStatus(std::vector<int>& trace, bool vote){

    int finalState = parse(trace, rootState);

    if(finalState == -1){
        if(vote){
            int lastKnown = partParse(trace, rootState);
            return merger.getNFutureAcc(lastKnown) > merger.getNFutureRej(lastKnown) ? 1 : 0;
        }else{
            return -1;
        }
    }else{
        return getAcceptanceStatus(finalState);
    }
}


/// @brief Checks for a given trace wether it is predicted correctly by the model
/// @param trace the trace which is to be predicted. Contains label information
/// @param style Encoding style of the trace. 1 -- is scientific standard
/// @param vote Wether a voting scheme is employed if the status of the trace is unkown. Will lable those as false otheriwse.
/// @return 0 -- for true negative; 1 -- for true positive; 2 -- for false negative; 3 -- for false positive
int InterfixGraph::evaluateAcc(const std::string& trace, int style, bool vote){

    std::vector<int> seq;
    auto lenacc = traceToSeq(trace, seq, style);
    
    int acc = getAcceptanceStatus(seq, vote);

    if(!lenacc.second){
        if(acc == 0 || acc == -1){
            return 0;
        }else{
            return 3;
        }
    }else{
        if(acc == 1){
            return 1;
        }else{
            return 2;
        }
    }
    
}

/// @brief Predict all traces in a file.
/// @param filepath Path to the file
/// @param details true negative count, true positive count, false negative count, false positive count, number of traces will be recorded here
/// @param style Encoding style of the file of traces. 1 is scientific standard
/// @param vote Wether voting scheme is employed for traces where the model is unsure about solution
/// @return Abolute accuracy on the given file
float InterfixGraph::testAccFile(const std::string& filepath, std::vector<int>& details, int style, bool vote){

    details = std::vector<int>(5, 0);

    std::ifstream fs(filepath);

    std::string line;

    if(std::getline(fs, line)){

        std::istringstream ss(line);
        int nTraces, nSymbols;
        ss >> nTraces >> nSymbols;

        details[4] = nTraces;

        int truePos = 0;
        int trueNeg = 0;
        int falsePos = 0;
        int falseNeg = 0;

        for(int i = 0; i < nTraces; i++){

            if(std::getline(fs, line)){

                int result = evaluateAcc(line, style, vote);

                switch (result)
                {
                case 0:
                    trueNeg++;
                    break;

                case 1:
                    truePos++;
                    break;

                case 2:
                    falseNeg++;
                    break;

                case 3:
                    falsePos++;
                    break;
                
                default:
                    break;
                }

            }else{
                std::cout << "Error: couldn't read line " << i << " from file " << filepath << "\n";
                break;
            }
        }

        details[0] = trueNeg;
        details[1] = truePos;
        details[2] = falseNeg;
        details[3] = falsePos;

        float acc = ((float) trueNeg + (float) truePos) / (float) nTraces;

        return acc;

    }else{
        std::cout << "Error: couldn't read from file " << filepath << "\n";
        return 0.0;
    }

}

/// @brief Create a node label for the .dot file for a given node. 
/// @param node 
/// @param redOnly If red only mode is set, red states will be coloured red, blue states blue, and rejecting states will be diamonds
/// @return Double circle for accepting state, circle with bars for rejecting state, regular circle for unknown states
std::string InterfixGraph::dotNodeType(const int node, bool redOnly){

    std::string modifier = "";

    modifier = " [";
    if(redOnly){
        if(getSColour(node) == 1){
            modifier = modifier + "color = red, style=filled";
        }else{
            modifier = modifier + "color = blue, style=filled";
        }
    }

    switch (getAcceptanceStatus(node)){

    case 1:
        if(redOnly){
            modifier = modifier + ", shape=doublecircle]";
            return modifier;
        }
        return " [shape = doublecircle]";

    case 0:
        if(redOnly){
            modifier = modifier + ", shape=diamond]";
            return modifier;
        }
        return " [shape = Mcircle]";
    
    default:
        if(redOnly){
            modifier = modifier + "]";
            return modifier;
        }
        return "";

    }
}

/// @brief Recursively explore (BFS) the graph and write it to a .dot file
/// @param os Out stream to which the graph is written
/// @param visited Vector if a node is already encountered (needs to be equal in lenght to number of states)
/// @param cNode the current node that is considered
/// @param name the name used for the current node
/// @param nameReg A vector mapping all nodes to their name if it is already known (length needs to be equal to number of states)
/// @param toExplore Queue of states that will be explored
/// @param right Whether the right cayley graph will be recorded in the .dot file
/// @param left Whether the left cayley graph will be recorded in the .dot file
/// @param redOnly Wheter only red states will be recorded in the .dot file
void InterfixGraph::recursiveDotWrite(std::ostream& os, std::vector<bool>& visited, int cNode, 
                                      std::string name, std::vector<std::string>& nameReg,
                                      std::queue<std::pair<int,std::string>>& toExplore,
                                      bool right, bool left, bool redOnly){
    
    std::string cNodeName;
    if(cNode == rootState){
        cNodeName = "root";
    }else{
        cNodeName = name;
    }

    int nextNode;

    if(right && (!redOnly || getSColour(cNode) == 1)){
        //Do BFS for shortest names possbile
        for(int i=0; i < nSymbols; i++){
            nextNode = followEdge(cNode, i);
            if(nextNode != -1){
                if(!visited[nextNode]){
                    std::string nextNodeName = name + std::to_string(i);
                    os << nextNodeName << dotNodeType(nextNode, redOnly) << "\n";
                    os << cNodeName << " -> " << nextNodeName << " [color = red, label = "<< i << "] \n";
                    toExplore.push(std::pair(nextNode, nextNodeName));
                    visited[nextNode] = true;
                    nameReg[nextNode] = nextNodeName;
                }else{
                    os << cNodeName << " -> " << nameReg[nextNode] << " [color = red, label = "<< i << "] \n";
                }
            }
        }
    }
    if(left && (!redOnly || getSColour(cNode) == 1)){
        for(int i=0; i < nSymbols; i++){
            nextNode = followEdge(cNode, mapSymb(i, true));
            if(nextNode != -1){
                if(!visited[nextNode]){
                    std::string nextNodeName = std::to_string(i) + name;
                    os << nextNodeName << dotNodeType(nextNode, redOnly) << "\n";
                    os << cNodeName << " -> " << nextNodeName << " [color = blue, label = "<< i << "] \n";
                    toExplore.push(std::pair(nextNode, nextNodeName));
                    visited[nextNode] = true;
                    nameReg[nextNode] = nextNodeName;
                }else{
                    os << cNodeName << " -> " << nameReg[nextNode] << " [color = blue, label = "<< i << "] \n";
                }
            }
        }
    }

}

/// @brief Write the graph to a .dot file
/// @param filepath Path to which the .dot file is to be written
/// @param right Whether the right Cayley graph will be included
/// @param left Wheter the left Cayley graph will be included
/// @param redOnly Whether only red states will be included
/// @return 1 on success; 0 on failure
int InterfixGraph::to_dot_file(const std::string& filepath, bool right, bool left, bool redOnly){

    std::ofstream fs(filepath);
    if(fs.is_open()){

        std::vector<bool> drawn(getNNodes(), false);
        std::string cLabel = "";

        fs << "digraph{ \n";

        fs << "root " << dotNodeType(rootState, redOnly) << "\n";
        
        drawn[rootState] = true;
        std::vector<std::string> nameReg(getNNodes());
        nameReg[rootState] = "root";
        std::queue<std::pair<int,std::string>> toExplore;
        toExplore.push(std::pair(rootState, ""));

        while(!toExplore.empty()){
            recursiveDotWrite(fs, drawn, (toExplore.front()).first, 
                              (toExplore.front()).second, nameReg, toExplore, right, left, redOnly);
            toExplore.pop();
        }
        

        fs << "} \n";

        fs.close();

    }else{
        std::cout << "Couldn't write to file: " << filepath << std::endl;
        return 0;
    }

    return 1;

}

/// @brief Writes the graph to a .txt file. Compatible with the read constructor
/// @param filepath Path to the file where it will be written
/// @return 1 on success; 0 on failure
int InterfixGraph::to_txt_file(const std::string& filepath){

    std::ofstream os(filepath);
    if(os.is_open()){

        std::queue<int> toExplore;
        std::vector<int> idMap(getNNodes(), -1);
        std::vector<bool> visited(getNNodes(), false);

        int newID = 0;

        os << merger.getNSets() << " " << nSymbols << "\n";
        toExplore.push(rootState);
        idMap[rootState] = 0;
        newID++;
        visited[rootState] = true;

        while(!toExplore.empty()){
            int cState = toExplore.front();
            int cName = idMap[cState];
            os << cName << " " << merger.getAcc(cState) << " " << merger.getFreq(cState) << "\n";
            for(int c=0; c < getNColours(); c++){
                int next = followEdge(cState, c);
                int nName = -1;
                if(next != -1){
                    if(!visited[next]){
                        visited[next] = true;
                        toExplore.push(next);
                        nName = newID;
                        idMap[next] = nName;
                        newID++;
                    }else{
                        nName = idMap[next];
                    }
                }
                os << cName << " " << c << " " << nName << " " << merger.getOutFreq(cState, c) << "\n";
            }
            toExplore.pop();
        }

    }else{
        std::cout << "Couldn't write to file: " << filepath << std::endl;
        return 0;
    }

    os.close();
    return 1;

}

