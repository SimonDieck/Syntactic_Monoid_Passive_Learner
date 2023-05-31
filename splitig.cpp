#include "splitig.h"
#include "interfixgraph.h"
#include "automaton.h"

#include <vector>
#include <string>

#include <algorithm>


SplitIG::SplitIG(int nSymbols, const int nStates) : InterfixGraph(){
    delayedInitalize(nSymbols, nStates);
}

SplitIG::SplitIG(const std::string& filepath) : InterfixGraph(filepath){
    reverseSymb = nSymbols;
}

SplitIG::SplitIG() : InterfixGraph(){
    reverseSymb = -2;
}

void SplitIG::delayedInitalize(int nSymbols, const int nStates){
    this->nSymbols = nSymbols + 1;
    rootState = 0;
    merger = UF(nStates, nSymbols + 1);
    //we will use the last open symbol id as the id for the reverse symbol
    reverseSymb = nSymbols; 
    ECGraph::delayedInitalize(nSymbols+1, nStates);
}

/// @brief Change the colour of a given state to red. Updates list of red states and the colour of all successors in the right graph.
/// @param state Id of the state that will be turned into a red state
void SplitIG::addRedState(int state){
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
    //Ignore the reverse symbol for colouring successors blue
    for(int c = 0; c < getNColours()-1; c++){
        int next = followEdge(i, c);
        if(next != -1){
            if(colour[next] == -1){
                colour[next] = 0;
            }
        }
    }
}

/// @brief Identify the blue state with the highest frequency
/// @return The id of the blue state with the highest frequency
int SplitIG::getMaxBlue(){
    bool hasFreq = merger.hasFreq();
    int i = -1;
    int max = -1;
    for(auto it = red.begin(); it != red.end(); ++it){
        i = merger.find(*it);
        for(int c = 0; c < getNColours()-1; c++){
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

bool SplitIG::checkTermination(){
    return merger.getNSets()/2 > red.size();
}

void SplitIG::enrichWithTrace(const std::string& trace, int style, bool interfix){
    finalised = false;
    std::vector<int> seq;
    auto lenAcc = traceToSeq(trace, seq, style);
    int length = lenAcc.first;
    bool accepted = lenAcc.second;

    std::vector<int> revSeq = std::vector<int>(seq);
    revSeq.push_back(reverseSymb);
    std::reverse(revSeq.begin(), revSeq.end());
    enrichWithSuffixTrace(revSeq, 0, rootState);

    for(int i=0; i < length; i++){
        enrichWithSuffixTrace(seq, i, rootState);
        std::vector<int> rightSeq = std::vector<int>(seq.begin()+i, seq.end());
        std::vector<int> corrLeftSeq = std::vector<int>(revSeq.begin(), revSeq.end()-i);
        int finalRightState = ECGraph::parse(rightSeq, rootState);
        int corrLeftState = ECGraph::parse(corrLeftSeq, rootState);
        ECGraph::adjustEdge(finalRightState, reverseSymb, corrLeftState);
    }

    //Make the correct left state an accepting state
    int finalState = ECGraph::parse(revSeq, rootState);
    accepted? acceptedStateList.push_back(finalState) : rejectedStateList.push_back(finalState);

}

int SplitIG::getAcceptanceStatus(std::vector<int>& trace, bool vote){
    if(!vote){
        std::vector<int> extendedTrace = std::vector<int>(trace);
        extendedTrace.push_back(reverseSymb);
        return InterfixGraph::getAcceptanceStatus(extendedTrace, vote);
    }else{
        std::vector<int> results(3, 0);
        int l = trace.size();
        for(int i=0; i <= l; i++){
            std::vector<int> rightTrace = std::vector<int>(trace.begin()+i, trace.end());
            std::vector<int> leftTrace = std::vector<int>(trace.begin(), trace.begin()+i);
            leftTrace.push_back(reverseSymb);
            std::reverse(leftTrace.begin(), leftTrace.end());
            rightTrace.insert(rightTrace.end(), leftTrace.begin(), leftTrace.end());
            int pred = InterfixGraph::getAcceptanceStatus(rightTrace, false);
            if(pred == -1){
                results[2]++;
            }else{
                results[pred]++;
            }
        }
        if(results[1] > results[0]){
            return 1;
        }else if(results[0] == 0){
            return -1;
        }else{
            return 0;
        }
    }
}