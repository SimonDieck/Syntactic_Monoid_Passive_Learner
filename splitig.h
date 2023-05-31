#ifndef SPLIT_INTERFIX_GRAPH
#define SPLIT_INTERFIX_GRAPH

#include "interfixgraph.h"

#include <vector>
#include <string>

class SplitIG: public InterfixGraph{

    int reverseSymb = -2;

    //float scoreMerge();

    void delayedInitalize(int nSymbols, const int nStates=1);

    void addRedState(int state);

    int getMaxBlue();

    bool checkTermination();

    public:

        SplitIG(int nSymbols, const int nStates=1);

        SplitIG(const std::string& filepath);

        SplitIG();

        void enrichWithTrace(const std::string& trace, int style=0, bool interfix=true);

        int getAcceptanceStatus(std::vector<int>& trace, bool vote=false);
};

#endif