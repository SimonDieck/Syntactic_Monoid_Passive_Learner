#ifndef UNION_FIND_CLASS
#define UNION_FIND_CLASS

#include <vector>
#include <stack>
#include <utility>

/// @brief A union-find data structure that works with edge-coloured graphs
class UF{

    /// @brief For each node the id of its current representative
    std::vector<int> reps;
    /// @brief For each node how deep the representative stack currenty is
    std::vector<int> depth;
    /// @brief For each node, the id of another node, that contributed an edge of a certain colour during merging, if it exists.
    std::vector<std::vector<int>> edgeMerge;
    /// @brief All pairs of nodes that were merged since the last time merges where accepted
    std::stack<std::pair<int,int>> mergeHist;
    /// @brief A list of all pairs that contributed new edges since the last time merges where accepted
    std::stack<std::pair<int,int>> edgeMHist;
    /// @brief For each node, whether it is currently accepting
    std::vector<int> acceptanceStatus;
    /// @brief For each node, how many of its futures are accpeting/rejecting (disregarding loops)
    std::vector<std::pair<int,int>> accFuture;
    /// @brief A history of all instances, where the acceptance status was updated due to merging
    std::stack<std::pair<int,int>> accHist;
    /// @brief For each node and all its outgoing edges their current frequency
    std::vector<std::pair<float,std::vector<float>>> frequencies;

    /// @brief How many sets there are currently in the Union Find structure
    int nSets;
    /// @brief How many colours edges can have
    int nColours=0;
    /// @brief Whether path shortening is employed for merging
    bool shorten=true;
    /// @brief Whether merges are always forced to happen in certain order regardless of depth
    bool forceleft=false;
    /// @brief Whether frequency information is tracked and updated during merging
    bool tracksFreq=false;

    public:
        UF();

        UF(int n);

        UF(int n, int c);

        UF(UF& uf);

        void addSet();

        int getNSets();

        int getMergedEdge(int state, int colour);

        void getAllReps(std::vector<int>& repList);

        bool hasFreq();

        void finaliseAcc(std::vector<int>& acc, std::vector<int>& rej);

        void initaliseAcc(int nNodes);

        void initaliseFreq(int i, int c);

        int getAcc(int i);

        void setAcc(int i, int status, bool track=true);

        int getNFutureAcc(int i);

        int getNFutureRej(int i);

        void incFutureAcc(int i, bool acc);

        float getFreq(int i);

        float getOutFreq(int i, int c);

        void incFreq(int i, float a=1.0);

        void incOutFreq(int i, int c, float a=1.0);

        int find(int p);

        void merge(int i, int j);

        void mergeEdge(int i, int c, int j);

        void mergeFreq(int i, int j);

        void undoMergeHist();

        void acceptMergeHist();

};

#endif