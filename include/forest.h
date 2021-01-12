#ifndef FOREST_H_
#define FOREST_H_

#include "tree.h"
#include "data.h"
#include <iostream>
#include "hyperparameters.h"
#include <string>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <libxml/tree.h>
#include <libxml/parser.h>
using namespace boost::numeric::ublas;

#ifdef USE_CUDA
#include "cudatemplates/devicememorypitched.hpp"
#include "cudatemplates/devicememorylinear.hpp"
#include "cudatemplates/hostmemoryheap.hpp"
#include "cuda/icg_rf_interface.cuh"
#endif

#ifdef WIN32
#pragma warning( disable : 4290 )
//warning C4290:  C++ exception specification ignored except to indicate a function
//                is not __declspec(nothrow)
#endif

class Forest
{
public:
    Forest(const HyperParameters &hp);
    Forest(const HyperParameters &hp, const std::string& forestFilename,bool recursive=true);
    Forest();
    ~Forest();

    xmlNodePtr save(FILE *fp) const;
 void evaConf(const matrix<float>& confidence,int m,int n,int width,int highth,int winm,int winn,int classnum);////////////////
    void train(const matrix<float>& data, const std::vector<int>& labels, bool use_gpu = false);
    void train(const matrix<float>& data, const std::vector<int>& labels, const std::vector<double>& weights, bool use_gpu = false);

    //matrix<float>
void eval(const matrix<float>& data, const std::vector<int>& labels, bool use_gpu = false);

    void trainAndTest1vsAll(const matrix<float>& dataTr, const std::vector<int>& labelsTr, const matrix<float>& dataTs, const std::vector<int>& labelsTs);

    std::vector<int> getPredictions() const { return m_predictions; };
    matrix<float> getConfidences() const { return m_confidences; };

    void save(const std::string &name = "default");
    void load(const std::string &name = "default");

    double oobe() const;

    std::vector<std::vector<int> > getPath(const matrix<float>& data, const int treeIndex) { return m_trees[treeIndex].getPath(data); }

    matrix<float> m_confidences;		// <-- moved from protected to public by ZEESHAN
    std::vector<int> m_predictions;

#ifdef USE_CUDA
    void createForestMatrix();
#endif

    void evalByCPU(const matrix<float>& data);	// <-- added by ZEESHAN

 protected:
    std::vector<Tree> m_trees;

#ifdef USE_CUDA
    Cuda::Array<float,2> *m_forest_d;
#endif

    HyperParameters m_hp;
   
    void writeError(const std::string& dataFileName, double error);

    void initialize(const int numSamples);

    void trainByCPU(const matrix<float>& data, const std::vector<int>& labels);
    void trainByCPU(const matrix<float>& data, const std::vector<int>& labels, const std::vector<double>& weights);
    void evalByCPU(const matrix<float>& data, const std::vector<int>& labels);

    void trainByGPU(const matrix<float>& data, const std::vector<int>& labels, const std::vector<double>& weights);
    void evalByGPU(const matrix<float>& data, const std::vector<int>& labels);

    double computeError(const std::vector<int>& labels);
    double computeError(const std::vector<int>& labels, const matrix<float>& confidences,
                        const std::vector<int>& voteNum);

};

#endif /* FOREST_H_ */
