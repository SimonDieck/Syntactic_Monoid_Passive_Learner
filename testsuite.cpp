#include "automaton.h"
#include "interfixgraph.h"
#include "testsuite.h"
#include "splitig.h"

#include <vector>
#include <unordered_set>
#include <iterator>

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <functional>

#include <random>
#include <cmath>
#include <numeric>

/// @brief Generates a trace for a composite Automaton, consisting of a and b
/// @param a 
/// @param b 
/// @param length Length of the trace that is generated
/// @return The trace as a string with format: membership, length, trace
std::string generateCompositeTrace(Automaton a, Automaton b, int length){

    auto gen = getRandomGen();
    std::uniform_int_distribution<int> splitDist(0,length);
    auto rSplit = std::bind(splitDist, gen);
    int burn = rSplit();

    int la = length - rSplit();
    int lb = length - la;

    std::string aTrace = a.generateTrace(la, 1);
    std::string bTrace = b.generateTrace(lb, 1);

    bool aAcc, bAcc;
    std::istringstream ssa(aTrace);
    std::istringstream ssb(bTrace);
    ssa >> aAcc >> la;
    ssb >> bAcc >> lb;

    bool acc = aAcc && bAcc;

    std::ostringstream oss;
    oss << acc << " " << length << " ";
    for(int i=0; i < la; i++){
        int next;
        ssa >> next;
        oss << next << " ";
    }
    for(int i=0; i < lb; i++){
        int next;
        ssb >> next;
        next += a.getNColours();
        oss << next << " ";
    }

    oss << "\n";

    return oss.str();

}

/// @brief Generates a file of traces for an Automaton which is a composition of a and b
/// @param filelocation Path to the file, where the traces will be written
/// @param a 
/// @param b 
/// @param nTraces How many traces the file will contain
/// @param minLength The minimum length for a trace
/// @param maxLength The maximum length for a trace
void generateFileOfCompositeTraces(std::string filelocation, Automaton a, Automaton b, 
                                   int nTraces, int minLength, int maxLength){
    std::ofstream fs(filelocation);
    std::string trace;

    auto gen = getRandomGen();
    std::uniform_int_distribution<int> dist(minLength,maxLength);
    auto rLength = std::bind(dist, gen);

    fs << nTraces;
    fs << " " << a.getNColours() + b.getNColours() << "\n";

    for(int i=0; i<nTraces; i++){
        trace = generateCompositeTrace(a, b, rLength());
        fs << trace;
    }

    fs.close();

}

/// @brief Generate an entire dataset of traces from different composited Automaton pairs
/// @param folderpath path to the folder, where the files will be created
/// @param nFiles how many files will be created
/// @param filesizes For each size it will generate a new file, with that many traces for each instance
void generateCompositeDataSet(std::string& folderpath, int nFiles, std::vector<int> filesizes){

    std::string endModel = "_model.txt";
    std::string endTraces = "_traces.txt";
    std::string endTrain = "_train.txt";
    std::string endTest = "_test.txt";

    auto gen = getRandomGen();
    std::uniform_int_distribution<int> Nnodesdist(3,12);
    std::uniform_int_distribution<int> nCdist(2,5);
    auto rNnodes = std::bind(Nnodesdist, gen);
    auto rNC = std::bind(nCdist, gen);
    //remove the first random variable from generator, as it generates deterministically due to some c++ quirks
    int burn = rNnodes();
    burn = rNC();

    for(int i=1; i <= nFiles; i++){
        int anC = rNC();
        int bnC = rNC();
        int anNodes = rNnodes();
        int bnNodes = rNnodes();
        std::uniform_int_distribution<int> anAccdist(1, anNodes);
        std::uniform_int_distribution<int> bnAccdist(1, bnNodes);
        auto raNacc = std::bind(anAccdist, gen);
        auto rbNacc = std::bind(bnAccdist, gen);
        burn = raNacc();
        burn = rbNacc();
        std::uniform_int_distribution<int> aAccdist(0, anC-1);
        std::uniform_int_distribution<int> bAccdist(0, bnC-1);
        auto raAcc = std::bind(aAccdist, gen);
        auto rbAcc = std::bind(bAccdist, gen);
        burn = raAcc();
        burn = rbAcc();
        int naAcc = raNacc();
        int nbAcc = rbNacc();
        std::vector<int> aAccList(naAcc, -1);
        std::vector<int> bAccList(nbAcc, -1);
        for(int j=0; j < naAcc; j++){
            aAccList[j] = raAcc();
        }
        for(int j=0; j < nbAcc; j++){
            bAccList[j] = rbAcc();
        }
        Automaton a = Automaton(anNodes, anC, 0, aAccList);
        Automaton b = Automaton(bnNodes, bnC, 0, bAccList);
        Automaton composite = Automaton(a, b);

        std::string model_path = folderpath + std::to_string(i) + endModel;
        std::ofstream os(model_path);
        composite.print(os);
        os.close();

        for(auto it=filesizes.begin(); it != filesizes.end(); ++it){
            std::string traces_path = folderpath + std::to_string(i) + "_" + std::to_string(*it) + "_" + endTraces;
            std::string train_path = folderpath + std::to_string(i) + "_" + std::to_string(*it) + "_" + endTrain;
            std::string test_path = folderpath + std::to_string(i) + "_" + std::to_string(*it) + "_" + endTest;

            generateFileOfCompositeTraces(traces_path, a, b, *it, 4, 32);
            split_file_into_train_test(traces_path, train_path, test_path);
        }
        

    }

}

/// @brief Generate a file of traces for Automaton a at the specified filelocation
/// @param filelocation 
/// @param a 
/// @param nTraces How many traces the file contains
/// @param minLength Minimum length of each trace
/// @param maxLength Maximum length of each trace
/// @param style Encoding of the traces. 1 -- is scientific standard
void generateFileOfTraces(std::string filelocation, Automaton a, int nTraces, int minLength, int maxLength, int style){

    std::ofstream fs(filelocation);
    std::string trace;

    auto gen = getRandomGen();
    std::uniform_int_distribution<int> dist(minLength,maxLength);
    auto rLength = bind(dist, gen);

    if(style == 0){
        a.print(fs);
    }
    fs << nTraces;
    if(style == 0){
        fs << "\n";
    }else{
        fs << " " << a.getNColours() << "\n";
    }

    for(int i=0; i<nTraces; i++){
        trace = a.generateTrace(rLength(), style);
        fs << trace;
    }

    fs.close();

}

/// @brief Loads the groundtruth of probabilities for a file of traces (used for pautomac)
/// @param filepath path to the file, where the ground truth is stored
/// @param dist Vector into which the ground truth is stored
void load_ground_truth(const std::string& filepath, std::vector<float>& dist){

    std::ifstream fs(filepath);

    if(fs.is_open()){

        std::string line;
        std::getline(fs, line);
        std::istringstream ss(line);

        int nTraces;
        ss >> nTraces;

        dist = std::vector<float>(nTraces, 0.0);

        for(int i=0; i < nTraces; i++){
            if(std::getline(fs, line)){

                std::istringstream ss(line);

                ss >> dist[i];

            }else{
                std:: cout << "Error: Couldn't read line " << i << " from file " << filepath << "\n";
            }
        }

        float sum = std::accumulate(dist.begin(), dist.end(), 0.0f);
        for(auto it=dist.begin(); it != dist.end(); ++it){
            *it = *it / sum;
        }

    }else{
        std::cout << "Error: Couldn't write to file " << filepath << "\n";
    }

}

/// @brief Given a vector containing ground truth probabilities and one containing predictions returns the perplexity of the two distributions 
/// @param truth vector containing ground truths, normalised to sum to 1.
/// @param predict vector containing all predictions, normalised to sum to 1.
/// @return the perplexity of the two distributions
float perplexity(std::vector<float>& truth, std::vector<float>& predict){

    if(truth.size() != predict.size()){
        std::cout << "Error: Number of instances between distributions differs. \n";
        return -1.0;
    }

    float dif = 0.0;

    for(int i = 0; i < (int) truth.size(); i++){
        dif += truth[i] * std::log2(predict[i]);
    }

    dif = -1 * dif;

    return std::exp2(dif);

}

/// @brief Splits a file of traces into a train and test file with an 80/20 split on unique traces
/// @param input Path to file of traces that will be split
/// @param out_train Path to the training file
/// @param out_test Path to the test file
/// @return 1 -- on failure, 0 -- on success
int split_file_into_train_test(const std::string& input, const std::string& out_train, const std::string& out_test){

    std::ifstream ins(input);
    std::ofstream os_train(out_train);
    std::ofstream os_test(out_test);

    //insert empty lines which are later overwritten with number of traces
    os_train << "                \n";
    os_test << "                \n";

    std::string line;

    if(std::getline(ins, line)){

        std::istringstream ss(line);

        int nTraces;
        int nSymbols;

        ss >> nTraces >> nSymbols; 

        //First pass is used to identify how many unique traces there are in the set
        std::unordered_set<std::string> uniqueTraces;

        for(int i=0; i < nTraces; i++){

            if(std::getline(ins, line)){

                uniqueTraces.insert(line);

            }else{
                std::cout << "Splitting stamina. First pass. \n";
                std::cout << "Error: couldn't read line " << i << " from file " << input << "\n";
                return 1;
            }

        }

        //split is based on unique traces. Duplicates of the same trace end in the same side of the split
        int nUnique = uniqueTraces.size();
        int nUTrain = (int) (0.8 * nUnique);
        auto eTrain = uniqueTraces.begin();
        std::advance(eTrain, nUTrain);
        std::unordered_set<std::string> train(uniqueTraces.begin(), eTrain);

        int nTrain = 0;
        int nTest = 0;

        //reset input stream and consume first line
        std::ifstream ins(input);
        std::getline(ins,line);

        for(int i=0; i < nTraces; i++){

            if(std::getline(ins, line)){

                if(train.count(line)){
                    os_train << line << "\n";
                    nTrain++;
                }else{
                    os_test << line << "\n";
                    nTest++;
                }

            }else{
                std::cout << "Splitting stamina. Second pass. \n";
                std::cout << "Error: couldn't read line " << i << " from file " << input << "\n"; 
                return 1;
            }

        }
        os_train.seekp(0);
        os_train << nTrain << " " << nSymbols;
        os_test.seekp(0);
        os_test << nTest << " " << nSymbols;

    }else{
        std::cout << "Error: couldn't read header from file " << input << " \n";
        return 1;
    }

    return 0;

}


/// @brief For a set, where a model was already trained this function loads the learned model and tests it again
/// @param dataModelPath Path to the learned models
/// @param nFiles How many files in the dataset will be retested
/// @param dataTestpath Path to the folder containing all the test files
/// @param modelName name of the model
/// @param outpath path to the folder, where results will be stored
/// @param vote TODO - voting doesn't work for retesting. Set to False
void retest_dataset(std::string& dataModelPath, int nFiles, std::string& dataTestpath, 
                    std::string& modelName, std::string& outpath, bool vote){
    
    for(int i=1; i<= nFiles; i++){

        std::string modelpath = dataModelPath + std::to_string(i) + modelName + "_model.txt";
        std::string inTest = dataTestpath + std::to_string(i) + "_test.txt";
        std::string outTest = outpath + std::to_string(i) + modelName + "_results.txt";

        retest_model(modelpath, inTest, outTest, vote);

    }


}

/// @brief Retests an already learned model
/// @param modelpath Path to the file containing the learned model
/// @param testfilepath path to the file containing all test traces
/// @param outpath path to the file where the results will be written
/// @param vote TODO- Voting doesn't work for retesting. Set to False
void retest_model(std::string& modelpath, std::string& testfilepath, std::string& outpath, bool vote){

    InterfixGraph model(modelpath);

    std::ofstream os(outpath);

    print_test(model, os, testfilepath, vote);

    os.close();

}

/// @brief Print the results for predicting a testfile with a given model
/// @param ifg the model that will predict the traces in the file
/// @param os stream to which results are written
/// @param testfilepath Path to the file containing the traces the model is tested on
/// @param vote Whether the voting scheme is employed for traces which lead to unknown transitions
/// @return The accuracy on the given testfile; -1.0 if an error occurs
float print_test(InterfixGraph& ifg, std::ofstream& os, std::string& testfilepath, bool vote, bool full){

    std::vector<int> details;

    float acc;

    acc = ifg.testAccFile(testfilepath, details, 1, vote);

    os << "Test accuracy: " << acc << "\n";
    os << "TN: " << details[0] << " TP: " << details[1] << 
          " FN: " << details[2] << " FP: " << details[3] << "\n";

    if(full){
        std::ifstream ts(testfilepath);

        std::string line;

        if(std::getline(ts, line)){

            std::istringstream ss(line);
            int nTraces, nSymbols;
            ss >> nTraces >> nSymbols;

            os << nTraces << " " << nSymbols << "\n";

            for(int j=0; j < nTraces; j++){

                if(std::getline(ts, line)){

                    int pred;
                    pred = ifg.evaluateAcc(line, 1, vote) % 2;

                    os << pred << "\n";

                }else{
                    std::cout << "Error: Couldn't read line " << j << " from file " << testfilepath << "\n";
                    ts.close();
                    return -1.0;
                }
            }

        }else{
            std::cout << "Error: Couldn't read test file " << testfilepath;
            ts.close();
            return -1.0;
        }

        ts.close();
    }
    
    return acc;

}

/// @brief Train models and test them for all files in a given dataset
/// @param in_path Path to the folder containing the dataset
/// @param out_path Path to the folder, where results will be stored
/// @param nFiles How many files of the dataset should be used
/// @param heuristic Which heuristic should be used to evaluate merges during state-merging
/// @param aut Wheter Automatons should be learned
/// @param sm Whether Syntactic Monoids should be learned
/// @param vote Wheter the voting scheme should be employed during evaluation
void train_test_dataset(std::string& in_path, std::string& out_path, int nFiles, int heuristic, 
                        bool aut, bool sm, bool vote, bool sig){
    
    std::string train = "_train.txt";
    std::string test = "_test.txt";

    std::string path_model_aut = "automaton\\";
    std::string path_model_sm = "syntactic_monoids\\";
    std::string path_model_sig = "split_graph\\";

    for(int i=1; i <= nFiles; i++){

        std::string in_train = in_path + std::to_string(i) + train;
        std::string in_test = in_path + std::to_string(i) + test;

        std::string out_aut = out_path + path_model_aut + std::to_string(i) + "_aut_model.txt";
        std::string out_sm = out_path + path_model_sm + std::to_string(i) + "_sm_model.txt";
        std::string out_sig = out_path + path_model_sig + std::to_string(i) + "_sig_model.txt";

        std::string out_test_aut = out_path + std::to_string(i)+"_aut_results.txt";
        std::string out_test_sm = out_path + std::to_string(i)+"_sm_results.txt";
        std::string out_test_sig = out_path + std::to_string(i)+"_sig_results.txt";

        float acc_aut = -1.0;
        float acc_sm = -1.0;
        float acc_sig = -1.0;

        if(aut){
            std::ofstream aos(out_test_aut);
            InterfixGraph pft;
            pft.enrichWithFileOfTraces(in_train, 1, true, false);
            pft.setScoreMode(heuristic);
            pft.collapse(2);
            pft.printInfo(std::cout);
            pft.printInfo(aos);
            pft.to_txt_file(out_aut);

            acc_aut = print_test(pft, aos, in_test, vote);
            aos.close();
        }

        if(sm){
            std::ofstream sos(out_test_sm);
            InterfixGraph ifg;
            ifg.enrichWithFileOfTraces(in_train, 1, true, true);
            ifg.setScoreMode(heuristic);
            ifg.collapse(2);
            ifg.printInfo(std::cout);
            ifg.printInfo(sos);
            ifg.to_txt_file(out_sm);

            acc_sm = print_test(ifg, sos, in_test, vote);
            sos.close();
        }

        if(sig){
            std::ofstream sios(out_test_sig);
            SplitIG sig;
            sig.enrichWithFileOfTraces(in_train, 1, true, false);
            sig.setScoreMode(heuristic);
            sig.collapse();
            sig.printInfo(std::cout);
            sig.printInfo(sios);
            sig.to_txt_file(out_sig);

            acc_sig = print_test(sig, sios, in_test, vote);
            sios.close();
        }
         
        std::cout << "Accuracy Automaton: " << acc_aut << "\n" << "Accuracy Syntactic Monoid: " << acc_sm <<"\n";
        std::cout << "Accuracy Split Graph: " << acc_sig << "\n";
        
    }

}

int main(){
    
    //generateCompositeDataSet(composite_inc_path, nIncComposite, filesizes);

    //generating the signature dataset
    /*
    auto gen1 = getRandomGen();
    std::uniform_int_distribution<int> nCdist(2,8);
    auto rnC = std::bind(nCdist, gen1);
    int burn = rnC();

    for(int i = 1; i <= 50; i++){

        int nC = rnC();
        int sigLBound = 18/nC;
        std::uniform_int_distribution<int> sigLdist(2, sigLBound);
        auto gen2 = getRandomGen();
        auto rSigL = std::bind(sigLdist, gen2);
        burn = rSigL();
        int sigLength = rSigL();
        std::uniform_int_distribution<int> colourDist(0, (nC-1));
        auto gen3 = getRandomGen();
        auto rC = std::bind(colourDist, gen3);
        burn = rC();
        std::vector<int> signature(sigLength, 0);
        for(int j=0; j < sigLength; j++){
            signature[j] = rC();
        }

        std::string sig_m_path = signature_path + std::to_string(i) + signature_model;
        std::string sig_traces_path = signature_path + std::to_string(i) + signature_traces;

        Automaton sigA(nC, signature, sigLength);
        std::ofstream os(sig_m_path);
        sigA.print(os);
        os.close();
        generateFileOfTraces(sig_traces_path, sigA, 5000, 1, nC*6, 1);

    }
    */

    std::string resultpath = ".\\test_results\\stamina\\stamina_agr\\";
    std::string smFolder = "syntactic_monoids\\";
    std::string autFolder = "automaton\\";
    std::string autName = "_aut";
    std::string smName = "_sm"; 
    std::string inpath = ".\\data\\staminadata\\";
    /*
    std::string modelpath = resultpath + autFolder;
    retest_dataset(modelpath, 25, inpath, autName, resultpath, false);
    modelpath = resultpath + smFolder;
    retest_dataset(modelpath, 25, inpath, smName, resultpath, false);
    */

    /*
    std::string input = inpath + "1_train.txt";
    InterfixGraph pft;
    pft.enrichWithFileOfTraces(input, 1, true, true);
    pft.setScoreMode(3);
    pft.collapse();
    */

    train_test_dataset(inpath, resultpath, 25, 2, false, false, true, true);

    return 0;
}