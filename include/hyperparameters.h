/*
 * hyperparametes.h
 *
 *  Created on: Jan 27, 2009
 *      Author: leisti
 */

#ifndef HYPERPARAMETES_H_
#define HYPERPARAMETES_H_
#include <string>

struct HyperParameters
{
  int numTrees;
  int maxTreeDepth;
  float bagRatio;
  float confThreshold;
  int numRandomFeatures;
  int numProjFeatures;
  int useRandProj;
  int useGPU;
  int useSubSamplingWithReplacement;
  int useSoftVoting;
  int useInfoGain;
  int numClasses;
  int verbose;
  int saveForest;
  int numLabeled;
  int numHistogramBins;
  int numTries;
int numTestNum;
int numTestNum1;
int numTestNum2;
int numTestNum3;
int weidu;
int width;
int highth;
int width1;
int highth1;
int width2;
int highth2;
int width3;
int highth3;
int width4;
int highth4;
int width5;
int highth5;
int width6;
int highth6;
int width7;
int highth7;
int width8;
int highth8;
int width9;
int highth9;
int width10;
int highth10;

int winm;
int winn;
int picnum;
int picnum0;
int partnum;
///////////////////////
  std::string className;
std::string className1;
std::string className2;
std::string className3;
std::string className4;
 std::string className5;
 std::string className6;
 std::string className7;
 std::string className8;
 std::string className9;
 std::string className10;

  std::string saveName;
  std::string savePath;
  std::string loadName;
  std::string trainData;
  std::string trainLabels;
std::string trainData1;
  std::string trainLabels1;
std::string trainData2;
  std::string trainLabels2;
std::string trainData3;
  std::string trainLabels3;
  std::string testData;
std::string testData1;
std::string testData2;
std::string testData3;
  std::string testLabels;
std::string testLabels1;
std::string testLabels2;
std::string testLabels3;
  int isExtreme;
};

#endif /* HYPERPARAMETES_H_ */
