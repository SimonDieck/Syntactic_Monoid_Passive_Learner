#include "unionfind.h"

#include <vector>
#include <numeric>
#include <algorithm>

/// @brief Default constructor. Information needs to be set afterwards
UF::UF(){
    nSets = 0;
}

/// @brief Constructor specifying the initial number of sets. Can't undo merges, when constructed this way
/// @param n Initial number of sets
UF::UF(int n){
    reps = std::vector<int>(n);
    std::iota(reps.begin(), reps.end(), 0);
    depth = std::vector<int>(n, 1);
    nSets = n;
}

/// @brief Full constructor specifying the initial number of sets as well as how many colours are in the multigraph. Can undo merges
/// @param n Initial number of sets
/// @param c Number of colours used in the edge-coloured graph
UF::UF(int n, int c){
    reps = std::vector<int>(n);
    edgeMerge = std::vector<std::vector<int>> (n, std::vector<int>(c, -1));
    std::iota(reps.begin(), reps.end(), 0);
    depth = std::vector<int>(n, 1);
    nSets = n;
    nColours = c;
    shorten = false;
    forceleft = true;
}

/// @brief Copy constructor. Will only copy data, that isn't related to unmerging. Thus merges can't be undone
/// @param uf Structure which is copied
UF::UF(UF& uf){
    this->nSets = uf.nSets;
    this->reps = std::vector<int>(uf.reps);
    this->depth = std::vector<int>(uf.depth);
}

/// @brief Add a new set with itself as representative
void UF::addSet(){
    reps.push_back(reps.size());
    depth.push_back(1);
    edgeMerge.push_back(std::vector<int>(nColours, -1));
    nSets++;
}

/// @brief Getter for number of sets
/// @return How many unique representatives are currently in the Union Find structure
int UF::getNSets(){
    return nSets;
}

/// @brief Obtain a list of all representatives
/// @param repList Representatives will be written into this vector
void UF::getAllReps(std::vector<int>& repList){
    for(int i = 0; i < (int) reps.size(); i++){
        if(reps[i] == i){
            repList.push_back(i);
        }
    }
}

/// @brief Given a list of rejecting and accepting states, makes sure that these properties are merged as well
/// @param acc List of all accepting states/representatives
/// @param rej List of all rejecting states/representatives
void UF::finaliseAcc(std::vector<int>& acc, std::vector<int>& rej){

    int n = reps.size();

    acceptanceStatus = std::vector<int>(n, -1);

    for(auto it = acc.begin(); it != acc.end(); ++it){
        acceptanceStatus[*it] = 1;
    }
    for(auto it = rej.begin(); it != rej.end(); ++it){
        acceptanceStatus[*it] = 0;
    }
}

/// @brief Initialises the vector containing acceptance status. Needs to be filled still
/// @param nNodes How many nodes there are in the graph. TODO -- This can be removed, as it is known to the UF
void UF::initaliseAcc(int nNodes){
    acceptanceStatus = std::vector<int>(nNodes, -1);
}

/// @brief Inititalise the vectors containing frequency information
/// @param i Number of states
/// @param c Number of colours in the edge-coloured graph
void UF::initaliseFreq(int i, int c){

    frequencies = std::vector<std::pair<float,std::vector<float>>> 
                        (i, std::pair<float,std::vector<float>>(0.0, std::vector<float>(c, 0.0)));
    accFuture = std::vector<std::pair<int,int>>(i, std::pair<int,int>(0,0));

    tracksFreq = true;
}

/// @brief Check whether the UF keeps track of frequency information
/// @return True if frequency is tracked. False otherwise
bool UF::hasFreq(){
    return tracksFreq;
}

/// @brief Get the acceptance status of a states representative
/// @param i state that is checked
/// @return 1 -- accepting; 0 -- rejecting; -1 -- Status unknown
int UF::getAcc(int i){
    return acceptanceStatus[find(i)];
}

/// @brief Set the acceptance status of a node to a specified value
/// @param i Node whose accpetance status is updated
/// @param status New acceptance status
/// @param track If true this change will be reverted if undomerge is called. If false the change is permanent
void UF::setAcc(int i, int status, bool track){
    if(track){
        accHist.push(std::pair<int,int>(i, acceptanceStatus[i]));
    }
    acceptanceStatus[i] = status;
}

/// @brief Check how many successor states of i (not considering loops) are accepting
/// @param i 
/// @return The number of accepting successors
int UF::getNFutureAcc(int i){
    return tracksFreq ? accFuture[i].first : 0;
}

/// @brief Check how many successor states of i (not considering loops) are rejecting
/// @param i 
/// @return The number of rejecting successors
int UF::getNFutureRej(int i){
    return tracksFreq ? accFuture[i].second : 0;
}

/// @brief Increase the number of accepting/rejecting futures of a state i by 1
/// @param i 
/// @param acc True if accepting futures should be increased, False, if rejecting futures should be increased
void UF::incFutureAcc(int i, bool acc){
    acc ? accFuture[i].first++ : accFuture[i].second++;
}

/// @brief Get the frequency of a state i
/// @param i id of the state
/// @return 
float UF::getFreq(int i){
    return tracksFreq ? frequencies[i].first : 0.0;
}

/// @brief Get the frequency of the edge from i of colour c
/// @param i id of the state
/// @param c id of the colour
/// @return 
float UF::getOutFreq(int i, int c){
    return tracksFreq ? frequencies[i].second[c] : 0.0;
}

/// @brief Increase the frequency of state i by a
/// @param i id of the state
/// @param a amount by which frequency is increased
void UF::incFreq(int i, float a){
    frequencies[i].first += a;
}

/// @brief Increase the frequency of an edge from i of colour c by a
/// @param i id of the state
/// @param c id of the colour
/// @param a amount by which the frequency is increased
void UF::incOutFreq(int i, int c, float a){
    frequencies[i].second[c] += a;
}

/// @brief If a merge introduced a new edge of colour to a state, this function will return it
/// @param state Id of the state which might have a new edge
/// @param colour Colour the edge would have
/// @return -1 if such an edge doesn't exist; The id of the node that contributed the new edge
int UF::getMergedEdge(int state, int colour){
    if(edgeMerge.empty()){
        return -1;
    }
    return edgeMerge[state][colour];
}

/// @brief Return the representative of a node p
/// @param p id of the node 
/// @return id of p's representative
int UF::find(int p){
    if(reps[p] == p){
        return p;
    }
    int root = find(reps[p]);
    if(shorten){
        reps[p] = root;
    }
    return root;
}

/// @brief Merges two nodes and all related information
/// @param i If forceleft is set, i will be the new representative. Otherwise depth is considered
/// @param j Id of the second state to be merged
void UF::merge(int i, int j){
    int repi = find(i);
    int repj = find(j);
    if(repi == repj){
        return;
    }
    if(depth[repj] < depth[repi] || forceleft){
        reps[repj] = repi;
        mergeHist.push(std::pair<int,int>(repj,repi));
        mergeFreq(repi, repj);
        depth[repi] = std::max(depth[repi], depth[repi]);
    }else if(depth[repj] > depth[repi]){
        reps[repi] = repj;
        mergeHist.push(std::pair<int,int>(repi,repj));
        mergeFreq(repj, repi);
    }else{
        reps[repj] = repi;
        mergeHist.push(std::pair<int,int>(repj,repi));
        mergeFreq(repi, repj);
        depth[repi]++;
    }
    nSets--;
    
}

/// @brief Makes j the origin for edges of colour c from i. Effectively adding the edge of colour c from j to i
/// @param i id of the first state
/// @param c id of the colour
/// @param j id of the second state
void UF::mergeEdge(int i, int c, int j){
    edgeMerge[i][c] = j;
    edgeMHist.push(std::pair<int,int>(i,c));
}

/// @brief Merges the frequencies of states. Frequency is accumulated into the record of i. So i should be representative of j afterwards
/// @param i id of i
/// @param j id of j
void UF::mergeFreq(int i, int j){
    if(!tracksFreq){
        return;
    }
    accFuture[i].first += accFuture[j].first;
    accFuture[i].second += accFuture[j].second;
    frequencies[i].first += frequencies[j].first;
    for(int c=0; c < nColours; c++){
        frequencies[i].second[c] += frequencies[j].second[c];
    }
}

/// @brief Reverts all changes to the structure made after the last time acceptMergeHist() was called
void UF::undoMergeHist(){
    int rep;
    int repee;
    while(!mergeHist.empty()){
        rep = mergeHist.top().first;
        repee = mergeHist.top().second;
        reps[rep] = rep;
        if(tracksFreq){
            accFuture[repee].first -= accFuture[rep].first;
            accFuture[repee].second -= accFuture[rep].second;
            frequencies[repee].first -= frequencies[rep].first;
            for(int c=0; c < nColours; c++){
                frequencies[mergeHist.top().second].second[c] -= frequencies[rep].second[c];
            }
        }
        mergeHist.pop();
        nSets++;
    }
    std::pair<int,int> edge;
    while(!edgeMHist.empty()){
        edge = edgeMHist.top();
        edgeMerge[edge.first][edge.second] = -1;
        edgeMHist.pop();
    }
    std::pair<int,int> acc;
    while(!accHist.empty()){
        acc = accHist.top();
        acceptanceStatus[acc.first] = acc.second;
        accHist.pop();
    }
}

/// @brief Makes all currently recorded changes irreversible. Undo merge will only undo until this point
void UF::acceptMergeHist(){
    shorten = true;
    while(!edgeMHist.empty()){
        edgeMHist.pop();
    }
    while(!accHist.empty()){
        accHist.pop();
    }
    while(!mergeHist.empty()){
        find(mergeHist.top().first);
        mergeHist.pop();
    }
    shorten = false;
}