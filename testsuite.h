#ifndef SMONOID_TESTSUITE
#define SMONOID_TESTSUITE

#include "automaton.h"
#include "interfixgraph.h"

#include <vector>

#include <string>
#include <iostream>

#include <random>

void load_ground_truth(const std::string& filepath, std::vector<float>& dist);

void retest_model(std::string& modelpath, std::string& testfilepath, std::string& outpath, bool vote=false);

void retest_dataset(std::string& dataModelPath, int nFiles, std::string& dataTestpath,
                    std::string& modelName, std::string& outpath, bool vote=false);

float print_test(InterfixGraph& ifg, std::ofstream& os, std::string& testfilepath, bool vote=false, bool full=true);

void train_test_dataset(std::string& in_path, std::string& out_path, int nFiles, int heuristic, 
                        bool aut, bool sm, bool vote=false, bool sig=false);

int split_file_into_train_test(const std::string& input, const std::string& out_train, const std::string& out_test);

float perplexity(std::vector<float>& truth, std::vector<float>& predict);

void generateCompositeDataSet(std::string& folderpath, int nFiles, std::vector<int> filesizes);

std::string generateCompositeTrace(Automaton a, Automaton b, int length);

void generateFileOfCompositeTraces(std::string filelocation, Automaton a, Automaton b, 
                                   int nTraces, int minLength, int maxLength);

void generateFileOfTraces(std::string filelocation, Automaton a, int nTraces, int minLength, int maxLength, int style=0);

#endif