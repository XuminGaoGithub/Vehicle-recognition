#include "/home/comin/libRandomForest_ZIA/include/forest.h"
//#include "forest.h"
#include <string.h>
#include<iostream>
#include <fstream>

#include <libxml/parser.h>
#include <boost/foreach.hpp>
#include<stdlib.h>
#include<stdio.h>
Forest::Forest(const HyperParameters &hp)
{
    m_hp = hp;
#ifdef USE_CUDA
    m_forest_d = NULL;
#endif

}

Forest::Forest(const HyperParameters &hp, const std::string& forestFilename,bool recursive) :

        m_hp( hp )
{cout<<"loading..."<<endl;////////////////////
	typedef boost::shared_ptr<NodeHyperPlane> PtrHyperPlane;
cout<<"10%"<<endl;////////////////////
	if (recursive == true)	{
cout<<"20%"<<endl;////////////////////
		xmlDocPtr forestDoc = xmlParseFile( forestFilename.c_str() );
cout<<"30%"<<endl;////////////////////
		if ( !forestDoc )
		{
		cout << "ERROR: no forest.xml file or wrong filename" << endl;
		return;
		}
cout<<"40%"<<endl;////////////////////
cout<<"50%"<<endl;////////////////////
		xmlNodePtr root = xmlDocGetRootElement( forestDoc );
cout<<"60%"<<endl;////////////////////
		if ( !root )
		{
		xmlFreeDoc( forestDoc );
		cout << "ERROR: no forest.xml file or wrong filename" << endl;
		return;
		}

		if ( xmlStrcmp( root->name, reinterpret_cast<const xmlChar*>( "randomforest" ) ) != 0 )
		{
		cout << "ERROR: no forest.xml file or wrong filename" << endl;
		cerr << "This doesn't seem to be classifier file..." << endl;
		xmlFreeDoc( forestDoc );
		return;
		}

		xmlNodePtr cur = root->xmlChildrenNode;
cout<<"70%"<<endl;////////////////////
		while ( cur != 0 )
		{
		if ( xmlStrcmp( cur->name, reinterpret_cast<const xmlChar*>( "tree" ) ) == 0 )
		{
		    m_trees.push_back( Tree( m_hp, cur ) );

		}
		cur = cur->next;
		}
cout<<"80%"<<endl;////////////////////
cout<<"90%"<<endl;////////////////////
		xmlFreeDoc( forestDoc );
cout<<"100%"<<endl;////////////////////
	}
	else {
cout<<"3"<<endl;////////////////////
   
		char str[500],str1[500],str2[500],str3[500],*pch,*pch1;
//		string str,str1,str2,str3;
		size_t pos;		
	        std::vector<int> bestFeatures;
     	        std::vector<float> bestWeights;
		std::vector<float> confidences;
	        float bestThreshold;
		bool isLeaf,leftchild;

		int tmpi,tmpi2,tmpi3,tmpi4,tmpi5,tmpi6,depth;
		float tmpf,tmpf2,tmpf3,tmpf4;
//		std::vector<PtrHyperPlane> stack_of_nodes;
		PtrHyperPlane curr_node;

		depth = 0;

		FILE *fp = fopen(forestFilename.c_str(),"r");
		fgets(str,450,fp);	// xml version 1.0
		fgets(str,450,fp);	// number of trees, already obtained from the configuration file
		printf("Number of trees = %d\n",hp.numTrees);

		for (int treeno=0;treeno!=hp.numTrees;treeno++)	{
			printf("Starting new tree...\n");
			fflush(stdout);
			fgets(str,450,fp);
			pch = strstr(str,"<tree ");	// Find <tree 
			sscanf(pch,"<tree num=\"%d\">\n",&tmpi);

			m_trees.push_back(Tree(m_hp));	// CREATING A TREE


			/////////////////////////////////////////
			// READING THE ROOT NODE
			/////////////////////////////////////////
			fgets(str,450,fp);		// node type
			fgets(str,450,fp);		// data

			bestFeatures.erase(bestFeatures.begin(),bestFeatures.end());
			bestWeights.erase(bestWeights.begin(),bestWeights.end());

			fgets(str,450,fp);		// bestfeatures
			for (int i=0;i!=hp.numProjFeatures;i++)	{
				fgets(str,450,fp);
				pch = strstr(str,"<feature value");	// Find <feature value
				sscanf(pch,"<feature value=\"%d\"/>\n",&tmpi);
							
				bestFeatures.push_back(tmpi);				
			}
			fgets(str,450,fp);		// /bestfeatures

			fgets(str,450,fp);		// bestweights
			for (int i=0;i!=hp.numProjFeatures;i++)	{
				fgets(str,450,fp);
				pch = strstr(str,"<weight value");	// Find <weight value
				sscanf(pch,"<weight value=\"%f\"/>\n",&tmpf);
							
				bestWeights.push_back(tmpf);				
			}
			fgets(str,450,fp);		// /bestweights
		
			fgets(str,450,fp);		// bestthreshold
			pch = strstr(str,"<bestthreshold");	// Find <bestthreshold
			sscanf(pch,"<bestthreshold value=\"%f\"/>\n",&bestThreshold);

			fgets(str,450,fp);		// </data>
			
			m_trees[treeno].m_rootNode = NodeHyperPlane::Ptr( new NodeHyperPlane(hp, 0, bestThreshold, bestFeatures, bestWeights) );	// root node
//			m_trees[treeno].m_rootNode.m_isLeaf = false;
			curr_node = boost::dynamic_pointer_cast<NodeHyperPlane>(m_trees[treeno].m_rootNode);
			curr_node->m_isLeaf = false;	
			curr_node->m_parent = boost::dynamic_pointer_cast<NodeHyperPlane>(curr_node);


			/////////////////////////////////////////
			// CHECK EACH NODE ONE BY ONE
			/////////////////////////////////////////
//cout<<"test"<<endl;
			while (true) {
				fgets(str,450,fp);		
				pch = strstr(str,"isLeaf=\"");	// NODE
				if (pch != NULL)	{
					pch1 = strstr(str,"true");
					if (pch1 != NULL)	isLeaf = true;
					else isLeaf=false;

					pch1 = strstr(str,"left"); 
					if (pch1 != NULL)	leftchild = true;
					else leftchild = false;

					if (isLeaf == false)	{				// NOT A LEAF
						fgets(str,450,fp);		// data

						bestFeatures.erase(bestFeatures.begin(),bestFeatures.end());
						bestWeights.erase(bestWeights.begin(),bestWeights.end());

						fgets(str,450,fp);		// bestfeatures
						for (int i=0;i!=hp.numProjFeatures;i++)	{
							fgets(str,450,fp);
							pch = strstr(str,"<feature value");	// Find <feature value
							sscanf(pch,"<feature value=\"%d\"/>\n",&tmpi);
							
							bestFeatures.push_back(tmpi);				
						}
						fgets(str,450,fp);		// /bestfeatures

						fgets(str,450,fp);		// bestweights
						for (int i=0;i!=hp.numProjFeatures;i++)	{
							fgets(str,450,fp);
							pch = strstr(str,"<weight value");	// Find <weight value
							sscanf(pch,"<weight value=\"%f\"/>\n",&tmpf);
							
							bestWeights.push_back(tmpf);				
						}
						fgets(str,450,fp);		// /bestweights
		
						fgets(str,450,fp);		// bestthreshold
						pch = strstr(str,"<bestthreshold");	// Find <bestthreshold
						sscanf(pch,"<bestthreshold value=\"%f\"/>\n",&bestThreshold);

						fgets(str,450,fp);		// </data>

						if (leftchild == true)	{
							depth++;
							curr_node->m_leftChildNode = NodeHyperPlane::Ptr( new NodeHyperPlane(hp, depth, bestThreshold, bestFeatures, bestWeights) );	// root node
							curr_node->m_leftChildNode->m_parent = boost::dynamic_pointer_cast<NodeHyperPlane>(curr_node);
							curr_node->m_isLeaf = false;	
							curr_node = boost::dynamic_pointer_cast<NodeHyperPlane>(curr_node->m_leftChildNode);
						}
						else {
							depth++;
							curr_node->m_rightChildNode = NodeHyperPlane::Ptr( new NodeHyperPlane(hp, depth, bestThreshold, bestFeatures, bestWeights) );	// root node
							curr_node->m_rightChildNode->m_parent = boost::dynamic_pointer_cast<NodeHyperPlane>(curr_node);
							curr_node->m_isLeaf = false;	
							curr_node = boost::dynamic_pointer_cast<NodeHyperPlane>(curr_node->m_rightChildNode);
						}



					}	
					else {
						pch1 = strstr(str,"label=");
						assert(pch1!=NULL);
						sscanf(pch1,"label=\"%d\"",&tmpi);

						if (leftchild == true)	{
							depth++;
							curr_node->m_leftChildNode = NodeHyperPlane::Ptr( new NodeHyperPlane(hp, depth) );	// root node
							curr_node->m_leftChildNode->m_parent = boost::dynamic_pointer_cast<NodeHyperPlane>(curr_node);
							curr_node->m_leftChildNode->m_nodeConf.erase(curr_node->m_leftChildNode->m_nodeConf.begin(),curr_node->m_leftChildNode->m_nodeConf.end());

							for (int i=0;i!=hp.numClasses;i++)	{
								fgets(str,450,fp);		// confidence
								pch1 = strstr(str,"conf=\"");
								sscanf(pch1,"conf=\"%f\"",&tmpf);
								curr_node->m_leftChildNode->m_nodeConf.push_back(tmpf);
								curr_node->m_leftChildNode->m_nodeLabel = tmpi;
								curr_node->m_leftChildNode->m_isLeaf = true;	
							}
						}
						else {
							depth++;
							curr_node->m_rightChildNode = NodeHyperPlane::Ptr( new NodeHyperPlane(hp, depth) );	// root node
							curr_node->m_rightChildNode->m_parent = boost::dynamic_pointer_cast<NodeHyperPlane>(curr_node);
							curr_node->m_rightChildNode->m_nodeConf.erase(curr_node->m_rightChildNode->m_nodeConf.begin(),curr_node->m_rightChildNode->m_nodeConf.end());

							for (int i=0;i!=hp.numClasses;i++)	{
								fgets(str,450,fp);		// confidence
								pch1 = strstr(str,"conf=\"");
								sscanf(pch1,"conf=\"%f\"",&tmpf);
								curr_node->m_rightChildNode->m_nodeConf.push_back(tmpf);
								curr_node->m_rightChildNode->m_nodeLabel = tmpi;
								curr_node->m_rightChildNode->m_isLeaf = true;	
							}
						}	
						fgets(str,450,fp);	// /node
					}	
				
				}
				else {
					pch = strstr(str,"</node>");	// /NODE
					if (pch != NULL)	{
						curr_node = boost::dynamic_pointer_cast<NodeHyperPlane>(curr_node->m_parent);
					}
					pch = strstr(str,"</tree>");	// /NODE
					if (pch != NULL)	{
						break;
					}
				}
			}
		}
		fclose(fp);
//cout<<"test"<<endl;
	}

}

void Forest::initialize(const int numSamples)
{

  m_confidences.resize(numSamples, m_hp.numClasses);
  m_predictions.resize(numSamples);

  for (int n = 0; n < numSamples; n++)
  {
    for (int m = 0; m < m_hp.numClasses; m++)
    {
        m_confidences(n, m) = 0;
    }
  }
}

xmlNodePtr Forest::save(FILE *fp) const
{
/*    xmlNodePtr node = xmlNewNode( NULL, reinterpret_cast<const xmlChar*>( "randomforest" ) );
    addIntProp( node, "numTrees", m_trees.size() );*/

	fprintf(fp,"<randomforest numTrees=\"%d\">\n",m_trees.size());

    //xmlAddChild( node, Configurator::conf()->saveConstants() );

    int i = 0;
    BOOST_FOREACH(Tree t, m_trees) {
      xmlNodePtr stageNode = t.save(fp,i);
/*      addIntProp( stageNode, "num", i );
      xmlAddChild( node, stageNode );*/
      i++;
    }

	fprintf(fp,"</randomforest>");

 //   return node;
}

void Forest::train(const matrix<float>& data, const std::vector<int>& labels, bool use_gpu)
{
    // Initialize
    initialize(data.size1());
//const << labels << endl;///////////////////////////////////////////////////////
    if (m_hp.useGPU || use_gpu){
        std::vector<double> weights(labels.size(),1);
        trainByGPU(data,labels, weights);
    }
    else {
        trainByCPU(data,labels);
    }
}

void Forest::train(const matrix<float>& data, const std::vector<int>& labels, const std::vector<double>& weights, bool use_gpu)
{
    // Initialize
    initialize(data.size1());

    if (m_hp.useGPU || use_gpu)
        trainByGPU(data,labels, weights);
    else
      trainByCPU(data,labels,weights);
}

Forest::~Forest()
{
    m_trees.clear();
#ifdef USE_CUDA
    if (m_forest_d != NULL)
        delete m_forest_d;
#endif
}
//////////////////////////////////step2eval/////////////////////////////////

//matrix<float> 
void Forest::eval(const matrix<float>& data, const std::vector<int>& labels, bool use_gpu)
{cout<<"step2eval:"<<endl;
  if (m_hp.useGPU || use_gpu)
    evalByGPU(data, labels);
  else
    evalByCPU(data, labels);
//cout << m_confidences(1, 1) << endl;///////////////////////////////
//return m_confidences;
}

void Forest::save(const std::string &name)
{
    std::string saveName;
    saveName = (name == "default") ? m_hp.saveName : name;

    // now save that stuff

}

void Forest::load(const std::string &name)
{
    std::string loadName;
    loadName = (name == "default") ? m_hp.loadName : name;

    // now load that stuff

}

#ifdef USE_CUDA
void Forest::createForestMatrix()
{
  // DEPRECATED !!
    // create a forest matrix for evaluation on the gpu
    int num_cols = 0;
    if (m_hp.useRandProj)
        num_cols = 3 + 2 * m_hp.numProjFeatures + 1 + m_hp.numClasses + 1;
    else
        num_cols = 3 + 2 + 1 + m_hp.numClasses + 1;

    int num_nodes_tree = (int)pow((float)2,m_hp.maxTreeDepth + 1) - 1;
    int num_rows = num_nodes_tree;

    float *forest_host = new float[num_cols * m_hp.numTrees * num_rows];

    if (m_hp.numTrees != m_trees.size())
        throw("The number of trees trained is different from the number configured");

    std::vector<Tree>::iterator tree_it = m_trees.begin();

    for (int tree_index = 0; tree_index < m_hp.numTrees; tree_index++) {
        matrix<float> tree_matrix(num_nodes_tree, num_cols, -1.0f);
        tree_it->getTreeAsMatrix(&tree_matrix, tree_index);

        // copy data
        for (int row = 0; row < num_rows; row++) {
            forest_host[row * m_hp.numTrees * num_cols + tree_index * num_cols] = (float) tree_index;
            for (int col = 1; col < num_cols; col++) {
                forest_host[row * m_hp.numTrees * num_cols + tree_index * num_cols + col] = tree_matrix(row, col);
            }
        }
        tree_it++;
    }

    if (m_hp.verbose)
        std::cout << "\ncreateForestMatrix(): forest matrix needs " << num_cols * m_hp.numTrees * num_rows * 4 / 1024 << " KB of memory\n";

    Cuda::HostMemoryReference<float, 2> forest_hmr(Cuda::Size<2>(num_cols * m_hp.numTrees, num_rows), forest_host);
    // DEBUG
    /*for (int row = 0; row < num_rows; row++) {
        for (int col = 0; col < num_cols * m_hp.numTrees; col++) {
            if (col%num_cols == 0)
                printf("\n");
            printf("%1.1f ", forest_host[row * num_cols * m_hp.numTrees + col]);
        }
        printf("\n");
    }*/

    //std::cout << "\n dimensions: " << forest_hmr.size[0] << "x" << forest_hmr.size[1] <<"\n";
    //std::cout << " max texture size: " << (int)pow(2.0f,16) << "x" << (int)pow(2.0f,15) <<"\n";
    if (forest_hmr.size[0] > pow(2.0f,14) || forest_hmr.size[1] > pow(2.0f,15)) {
        std::cout << "\nForest is larger than maximum cuda texture dimensions: \n";
        std::cout << "forest: w: " << forest_hmr.size[0] << ", h: " << forest_hmr.size[1];
        std::cout << " maximum texture size: w: " << pow(2.0f,14) << ", h: " << pow(2.0f,15);

    }
    if (m_forest_d != NULL)
        delete m_forest_d;
    m_forest_d = new Cuda::Array<float,2>(forest_hmr);

    delete[] forest_host;
}
#endif

// CPU penne code only below this line
void Forest::trainByCPU(const matrix<float>& data, const std::vector<int>& labels)
{
    std::vector<int> outOfBagVoteCount(m_hp.numLabeled, 0);
    matrix<float> outOfBagConfidences(m_hp.numLabeled,m_hp.numClasses);
    for (int i = 0; i < m_hp.numLabeled; i++)
    {
        for ( int j = 0; j < m_hp.numClasses; j++)
        {
            outOfBagConfidences(i,j) = 0.0;
        }
    }

    std::vector<int>::iterator preItr, preEnd;
    HyperParameters tmpHP = m_hp;
    tmpHP.verbose = false;

    if (m_hp.verbose)
    {
        cout << "Training a random forest with " << m_hp.numTrees << " , grab a coffee ... " << endl;
        cout << "\tTree #: ";
    }

    m_trees.clear();
    for (int i = 0; i < m_hp.numTrees; i++)
    {
        if (m_hp.verbose && !(10*i%m_hp.numTrees))
        {
            cout << 100*i/m_hp.numTrees << "% " << flush;
        }

        Tree t(tmpHP);
        t.train(data,labels, m_confidences, outOfBagConfidences, outOfBagVoteCount);
        m_trees.push_back(t);
    }

    if (m_hp.verbose)
    {
        cout << " Done!" << endl;
    }

    double error = computeError(labels);
    double outOfBagError = computeError(labels, outOfBagConfidences, outOfBagVoteCount);

    if (m_hp.verbose)
    {
        cout << "\tTraining error = " << error << ", out of bag error = " << outOfBagError << endl;
    }
}

void Forest::trainAndTest1vsAll(const matrix<float>& dataTr, const std::vector<int>& labelsTr,
                                const matrix<float>& dataTs, const std::vector<int>& labelsTs)
{
  int trueNumClasses = m_hp.numClasses;
  matrix<float> testProb(dataTs.size1(), trueNumClasses);
  m_hp.numClasses = 2;
  bool verbose = m_hp.verbose;
  m_hp.verbose = false;

  // 1 vs All loop
  if (verbose) {
    cout << "Training and testing with 1-vs-all for " << trueNumClasses << " classes ..." << endl;
  }

  std::vector<int> tmpLabelsTr(dataTr.size1()), tmpLabelsTs(dataTs.size1());
  matrix<float> tmpConf;
  double error;
  for (int nClass = 0; nClass < trueNumClasses; nClass++) {
    if (verbose) {
      cout << "\tClass #: " << nClass + 1 << "  ... ";
    }
    initialize(dataTr.size1());
    std::vector<double> tmpWeights(dataTr.size1(), 1.0);

    for (int nSamp = 0; nSamp < (int) dataTr.size1(); nSamp++) {
      if (labelsTr[nSamp] == nClass) {
        tmpLabelsTr[nSamp] = 1;
        tmpWeights[nSamp] = trueNumClasses;
      }
      else {
        tmpLabelsTr[nSamp] = 0;
      }
    }
    if (m_hp.useGPU)
        trainByGPU(dataTr, tmpLabelsTr, tmpWeights);
    else
        trainByCPU(dataTr, tmpLabelsTr, tmpWeights);


    for (int nSamp = 0; nSamp < (int) dataTs.size1(); nSamp++) {
      if (labelsTs[nSamp] == nClass) {
        tmpLabelsTs[nSamp] = 1;
      }
      else {
        tmpLabelsTs[nSamp] = 0;
      }
    }
    if (m_hp.useGPU)
      evalByGPU(dataTs, tmpLabelsTs);
    else
      evalByCPU(dataTs, tmpLabelsTs);

    error = computeError(tmpLabelsTs);
    if (verbose) {
      cout << " Classification error = " << error << endl;
    }

    for (int nSamp = 0; nSamp < (int) dataTs.size1(); nSamp++) {
      testProb(nSamp, nClass) = m_confidences(nSamp, 1);
    }
  }

  // Make the decision
  int bestClass;
  double bestConf, testError = 0;
  std::vector<double> preCount(trueNumClasses, 0), trueCount(trueNumClasses, 0);
  for (int nSamp = 0; nSamp < (int) dataTs.size1(); nSamp++) {
    bestClass = 0;
    bestConf = 0;
    for (int nClass = 0; nClass < trueNumClasses; nClass++) {
      if (bestConf < testProb(nSamp, nClass)) {
        bestConf = testProb(nSamp, nClass);
        bestClass = nClass;
      }
    }

    m_predictions[nSamp] = bestClass;

    if (labelsTs[nSamp] != bestClass) {
      testError++;
    }

    preCount[bestClass]++;
    trueCount[labelsTs[nSamp]]++;
  }

  for (int n = 0; n < trueNumClasses; n++) {
    preCount[n] /= (trueCount[n] + 1e-6);
  }

  if (verbose)
    {
      cout << "\tForest Test Error = " << testError/((double) dataTs.size1()) << endl;
      cout << "\tPrediction ratio: ";
      dispVector(preCount);
    }

  m_hp.numClasses = trueNumClasses;
}

void Forest::trainByCPU(const matrix<float>& data, const std::vector<int>& labels, const std::vector<double>& weights)
{
    std::vector<int> outOfBagVoteCount(data.size1(), 0);
    matrix<float> outOfBagConfidences(data.size1(),m_hp.numClasses);
    for ( unsigned int i = 0; i < data.size1(); i++)
    {
        for ( int j = 0; j < m_hp.numClasses; j++)
        {
            outOfBagConfidences(i,j) = 0.0;
        }
    }

    HyperParameters tmpHP = m_hp;
    tmpHP.verbose = false;

    if (m_hp.verbose)
    {
        cout << "Training a random forest with " << m_hp.numTrees << " trees, grab a coffee ... " << endl;
        cout << "\tTree #: ";
    }

    m_trees.clear();
    for (int i = 0; i < m_hp.numTrees; i++)
    {
        if (m_hp.verbose && !(10*i%m_hp.numTrees))
        {
            cout << 100*i/m_hp.numTrees << "% " << flush;
        }

        Tree t(tmpHP);
        t.train(data,labels, weights, m_confidences, outOfBagConfidences, outOfBagVoteCount);
        m_trees.push_back(t);
    }
    if (m_hp.verbose)
    {
        cout << " Done!" << endl;
    }

    double error = computeError(labels);
    double outOfBagError = computeError(labels, outOfBagConfidences, outOfBagVoteCount);

    if (m_hp.verbose)
    {
        cout << "\tTraining error = " << error << ", out of bag error = " << outOfBagError << endl;
    }
}


/////////////////////////////////step3evalByCPU/////////////////////////////////



void Forest::evalByCPU(const matrix<float>& data, const std::vector<int>& labels)
{cout<<"step3evalByCPU:"<<endl;
    // Initialize
    initialize(data.size1());//m_confidences=0
int numm=1;
    BOOST_FOREACH(Tree t, m_trees) {//class Tree     m_trees?
cout<<"treenum:"<<numm++<<endl;
      t.eval(data, labels, m_confidences);	// <-- EVALUATING THE TESTPOINTS WITH EACH INDIVIDUAL TREE
    }

    // divide confidences by number of trees
    m_confidences *= (1.0f / m_hp.numTrees);
//cout<<"m_hp.numTrees:"<<m_hp.numTrees<<endl<<"m_confidences:"<<m_confidences(1,1)<<endl;
    int bestClass;
    double bestConf;
    for (int nSamp = 0; nSamp < (int) data.size1(); nSamp++) {
      bestConf = 0;
      bestClass = 0;
      for(int nClass = 0; nClass < m_hp.numClasses; nClass++) {
        if (m_confidences(nSamp, nClass) > bestConf) {
          bestConf = m_confidences(nSamp, nClass);
	  bestClass = nClass;
        }
	//cout << m_confidences(nSamp, nClass) << " ";
      }
     // cout << endl;
/////////////////////////////////////////////////////////
      m_predictions[nSamp] = bestClass;
    }

    double error = computeError(labels);

    if (m_hp.verbose)
    {
        cout << "\tForest test error = " << error << endl;
    }
}

void Forest::evalByCPU(const matrix<float>& data)   // <-- ADDED BY ZEESHAN
{
    // Initialize
    initialize(data.size1());

    BOOST_FOREACH(Tree t, m_trees) {
      t.eval(data, m_confidences);	// <-- EVALUATING THE TESTPOINTS WITH EACH INDIVIDUAL TREE
    }
//cout<<"m_hp.numClasses:"<<m_hp.numClasses<<endl;
    // divide confidences by number of trees
    m_confidences *= (1.0f / m_hp.numTrees);

/*    int bestClass;
    double bestConf;
    for (int nSamp = 0; nSamp < (int) data.size1(); nSamp++) {
      bestConf = 0;
      bestClass = 0;
      for(int nClass = 0; nClass < m_hp.numClasses; nClass++) {
        if (m_confidences(nSamp, nClass) > bestConf) {
          bestConf = m_confidences(nSamp, nClass);
	  bestClass = nClass;
        }
//	cout << m_confidences(nSamp, nClass) << " ";
      }
//      cout << endl;

      m_predictions[nSamp] = bestClass;
    }

    double error = computeError(labels);

    if (m_hp.verbose)
    {
        cout << "\tForest test error = " << error << endl;
    }*/
}



double Forest::computeError(const std::vector<int>& labels)
{
    int bestClass, nSamp = 0;
    float bestConf;
    double error = 0;
cout<<"m_hp.numClasses:"<<m_hp.numClasses<<endl;
ofstream fout;
fout.open("tmp/confidence.txt");
  //cout<<labels[1]<<endl;///////////////////////////////////////////////////////
    BOOST_FOREACH(int pre, m_predictions) {
      bestClass = 0;
      bestConf = 0;
//if (nSamp==2)
//{
//cout<<"nSamp"<<nSamp<<endl;
//cout<<"label"<<labels[nSamp]<<endl;
//}
      for (int nClass = 0; nClass < (int) m_hp.numClasses; nClass++)

        {  if (nSamp==2)
{
//cout<<"test"<<m_confidences(nSamp, nClass)<<endl;///////////////////confidences
}
         fout<<m_confidences(nSamp, nClass)<<" ";
          if (m_confidences(nSamp, nClass) > bestConf)
            { 
              bestClass = nClass;
              bestConf = m_confidences(nSamp, nClass);
            }
        }
       
      pre = bestClass;
//cout<<pre<<endl;//////////////////////////////////////////////test label results
      if (bestClass != labels[nSamp])
        { //cout <<"1"<<endl;///////////////////////////////////////////////////////////////////////
          error++;
        }

      nSamp++;
fout<<endl;
    }
fout.close();
    error /= (double) m_predictions.size();

    return error;
}

double Forest::computeError(const std::vector<int>& labels, const matrix<float>& confidences,
                            const std::vector<int>& voteNum)
{
    double error = 0, bestConf;
    int sampleNum = 0, bestClass;
    for (int nSamp = 0; nSamp < (int) confidences.size1(); nSamp++)
    {
        bestClass = 0;
        bestConf = 0;
        if (voteNum[nSamp])
        {    
            sampleNum++;
            for (int nClass = 0; nClass < (int) m_hp.numClasses; nClass++)
            {
                if (confidences(nSamp, nClass) > bestConf)
                {
                    bestClass = nClass;
                    bestConf = confidences(nSamp, nClass);
                }
            }

            if (bestClass != labels[nSamp])
            {
                error++;
            }
        }
    }
    error /= (double) sampleNum;

    return error;
}


// GPU spaghetti code only below this line
void Forest::trainByGPU(const matrix<float>& data, const std::vector<int>& labels, const std::vector<double>& weights)
{
#ifdef USE_CUDA
    // test parameters
    if (!m_hp.useRandProj)
        m_hp.numProjFeatures = 1;

    // copy input data to device
    size_t cols = data.size2();
    size_t rows = data.size1();
    Cuda::HostMemoryReference<float, 2> data_hmr(Cuda::Size<2>(cols, rows), (float*)&data.data()[0]);
    Cuda::Array<float, 2> data_dmp(data_hmr);

    // copy labels to device
    Cuda::HostMemoryHeap<int, 1> labels_hmr(Cuda::Size<1>(labels.size()));
    Cuda::HostMemoryHeap<float, 1> weights_hmr(Cuda::Size<1>(labels.size()));
    for (int i = 0; i < labels.size(); i++){
        labels_hmr.getBuffer()[i] = labels[i];
        weights_hmr.getBuffer()[i] = weights[i];
    }
    Cuda::DeviceMemoryLinear<int, 1> labels_dml(labels_hmr);
    Cuda::DeviceMemoryLinear<float, 1> weights_dml(weights_hmr);

    // create forest structure
    int num_cols = iGetNumTreeCols(m_hp.numProjFeatures, m_hp.numClasses);
    Cuda::DeviceMemoryPitched<float,2>* forest_dmp;

    float oobe;
    iTrainForest(data_dmp, labels_dml, weights_dml, &forest_dmp, &oobe, data_dmp.size[1],
        m_hp.numClasses, m_hp.numTrees, m_hp.maxTreeDepth, num_cols, m_hp.numRandomFeatures, m_hp.numProjFeatures, m_hp.bagRatio);

    if (m_forest_d != NULL)
        delete m_forest_d;
    //printf("forest needs %3.3f KB of memory", forest_dmp->stride[0] * forest_dmp->size[1] * sizeof(float) / 1024.0f);
    printf(" oobe: %3.3f ",oobe);
    m_forest_d = new Cuda::Array<float,2>(*forest_dmp);
    delete forest_dmp; // FIXXME, one could directly pass the array
#endif
}

void Forest::evalByGPU(const matrix<float>& data, const std::vector<int>& labels)
{
#ifdef USE_CUDA
    if (!m_hp.useRandProj)
        m_hp.numProjFeatures = 1;

    // initialize
    initialize(data.size1());

    // copy input data to device
    size_t cols = data.size2();
    size_t rows = data.size1();
    Cuda::HostMemoryReference<float, 2> data_hmr(Cuda::Size<2>(cols, rows), (float*)&data.data()[0]);
    Cuda::Array<float, 2> data_dmp(Cuda::Size<2>(cols, rows));
    Cuda::copy(data_dmp, data_hmr);

    // create output structures
    size_t num_samples = rows;
    Cuda::DeviceMemoryPitched<float,2> confidences_dmp(Cuda::Size<2>(num_samples, m_hp.numClasses));
    Cuda::DeviceMemoryLinear1D<float> predictions_dmp(num_samples);

    bool soft_voting = false;
    if (m_hp.useSoftVoting) soft_voting = true;

    int num_tree_cols = iGetNumTreeCols(m_hp.numProjFeatures, m_hp.numClasses);

    iEvaluateForest(*m_forest_d, data_dmp, &confidences_dmp, &predictions_dmp,
    num_samples, m_hp.numTrees, m_hp.maxTreeDepth, num_tree_cols, m_hp.numClasses, m_hp.numProjFeatures, soft_voting);

    // copy output to host memory
    Cuda::HostMemoryHeap<float,2> confidences_hmh(confidences_dmp);
    Cuda::HostMemoryHeap<float,1> predictions_hmh(predictions_dmp);
    float* confidences_hmhp = confidences_hmh.getBuffer();
    float* predictions_hmhp = predictions_hmh.getBuffer();

    // assign confidences and predictions
    for (unsigned int nSamp = 0; nSamp < num_samples; nSamp++) {
        for (int nClass = 0; nClass < m_hp.numClasses; nClass++)
            m_confidences(nSamp,nClass) = confidences_hmhp[nClass * num_samples + nSamp];
        m_predictions[nSamp] = (int) predictions_hmhp[nSamp];
    }

    double error = computeError(labels);

    if (m_hp.verbose)
    {
        cout << "\tTest error = " << error << endl;
    }
#endif
}

void Forest::writeError(const std::string& fileName, double error)
{
    std::ofstream myfile;
    myfile.open(fileName.c_str(), std::ios::app);
    myfile << error << endl;
    myfile.flush();
    myfile.close();
}

void Forest::evaConf(const matrix<float>& confidence,int m,int n,int width,int highth,int winm,int winn,int classnum)
   { double bestconfidence[m];
     double bestclass[m];
//cout<<"1"<<endl;
for(int i=0;i<m;i++)
{
    //cout<<i<<endl;
     bestconfidence[i]=0;
     bestclass[i]=0;
    for(int j=0;j<n;j++)
      { 
        if(confidence(i,j)>bestconfidence[i])
        {    bestconfidence[i]=confidence(i,j);
            bestclass[i]=j;
        }
    }
}



//////////////////////nonzero
int num=0;
//std::vector<int> nonzeonum;

//std::vector<int> nonzeobestclass;
//std::vector<float> nonzeobestconfidence;
int nonzeonum[m];
int nonzeobestclass[m];
float nonzeobestconfidence[m];
matrix<float> nonzeroconfidence;
nonzeroconfidence.resize(m,n);
for(int ii=0;ii<m;ii++)
{
    if(bestclass[ii]!=0)
    {
        
        
        nonzeobestclass[num]=bestclass[ii];
        nonzeonum[num]=ii+1;
        nonzeobestconfidence[num]=bestconfidence[ii];
//cout<<nonzeobestclass[num]<<endl;
       for(int k=0;k<n;k++)
       {//nonzeroconfidence(0,0)=1;
 //cout<<confidence(0,k)<<endl;
        nonzeroconfidence(num,k)=confidence(ii,k);

       }
       num=num+1;

    }

}
//for(int k=0;k<num;k++)
//cout<<"nonzeonum"<<nonzeonum[k]<<endl;
//cout<<"3"<<endl;
/////////compute xylabel

float slide=8;
int windowl=width*2+1;
int windowh=highth*2+1;
int row,col;
//std::vector<int> Row;
int Row[num];
int Col[num];
float Rowlabel[num];
float Collabel[num];
//float temp[num];
//std::vector<int> Col;
//std::vector<int> Rowlabel;
//std::vector<int> Collabel;
row=floor(float(winm-windowh)/slide);
col=floor(float(winn-windowl)/slide);
//cout<<row<<endl<<col<<endl;
ofstream foutcol;
foutcol.open("tmp/col.txt");
ofstream foutrow;
foutrow.open("tmp/row.txt");
ofstream foutnonzero;
foutnonzero.open("tmp/nonzero.txt");
  for(int i=0;i<num;i++)
{//temp[i]=float(nonzeonum[i])/float(col);

Row[i]=floor(float(nonzeonum[i])/float(col))+1;
//cout<<temp[i]<<endl<<floor(nonzeonum[i]/col)<<endl;
 Col[i]=round((float(nonzeonum[i])/float(col)-floor(float(nonzeonum[i])/float(col)))*col);
//cout<<"Row"<<Row[i]<<endl;
//cout<<"Col"<<Col[i]<<endl;

foutnonzero<<nonzeonum[i]<<endl;
 }

 for(int i=0;i<num;i++)
{
     if(Col[i]==0)
{
         Row[i]=Row[i]-1;
         Col[i]=col;
    }
}

  for(int i=0;i<num;i++)
{

foutcol<<Col[i]<<endl;
foutrow<<Row[i]<<endl;
}
foutcol.close();
foutnonzero.close();foutrow.close();
//char *label,*classification;
//sprintf(label,"%d",classnum);
//sprintf(classification,"%d",classnum+100);
//cout<<"label:"<<label<<endl<<"classification:"<<classification<<endl;
char title[20]="label/label";
char classtitle[20]="label/class";
char confititle[25]="label/confidence";
char str[35];
char str1[10]=".txt";
int numtest=100;
//_itoa_s(numtest,str,10);

sprintf(str,"%d",classnum);
strcat(str,str1);
strcat(title,str);
strcat(classtitle,str);
strcat(confititle,str);
cout<<"str:"<<str<<endl;
//cout<<"confidence:"<<m_confidences(0,0)<<endl;

ofstream foutconf;
foutconf.open(confititle);
int nSamp=0;
 BOOST_FOREACH(int pre, m_predictions) {
for (int nClass = 0; nClass < (int) m_hp.numClasses; nClass++)

{  

        foutconf<<m_confidences(nSamp, nClass)<<" ";
        
 }
       
 
    nSamp++;
foutconf<<endl;
    }
foutconf.close();







ofstream foutlabel;
foutlabel.open(title);
//foutlabel.open(label);
ofstream foutclass;
foutclass.open(classtitle);
ofstream foutcollabel;
foutcollabel.open("tmp/collabel.txt");
ofstream foutrowlabel;
foutrowlabel.open("tmp/rowlabel.txt");
//foutclass.open(classification);
 for(int i=0;i<num;i++)
{Rowlabel[i]=(Row[i]-1)*slide+(windowh-1)/2;
Collabel[i]=(Col[i]-1)*slide+(windowl-1)/2;
//cout<<"Row"<<Row[i]<<endl;
//cout<<"Col"<<Col[i]<<endl;
//cout<<"Rowlabel"<<Rowlabel[i]<<endl;
//cout<<"Collabel"<<Collabel[i]<<endl;

foutcollabel<<Collabel[i]<<endl;
foutrowlabel<<Rowlabel[i]<<endl;

} 


foutcollabel.close();
foutrowlabel.close();


//std::vector<int> Cnonzerobestclass;
//std::vector<float> Cnonzeobestconfidence;
///std::vector<int> cRowlabel;
//std::vector<int> cCollabel;
int Cnonzerobestclass[num];
float Cnonzeobestconfidence[num];
float cRowlabel[num];
float cCollabel[num];
 float part1schp[4];
 float part2schp[4];
 float  part3schp[4];
 int partloop;
int global=1;
//choose
 int k=0;
cout<<"part1schp:"<<endl;

for(int i=0;i<4;i++)
{

cout<<part1schp[i]<<endl;
}

cout<<"part2schp:"<<endl;

for(int i=0;i<4;i++)
{

cout<<part2schp[i]<<endl;
}

cout<<"part3schp:"<<endl;

for(int i=0;i<4;i++)
{

cout<<part3schp[i]<<endl;
}

//cout<<"partloop:"<<partloop<<endl;


if (global==1)
{

 for(int i=0;i<num;i++)
 {  
//if(Rowlabel[i]<winm/1.5&&Rowlabel[i]>winm*0.33&&Collabel[i]>winn/3&&Collabel[i]<winn*1) 
 //{
         Cnonzerobestclass[k]=nonzeobestclass[i];
//cout<<"Cnonzerobestclass"<<Cnonzerobestclass[k]<<endl;
         cRowlabel[k]=Rowlabel[i];
         cCollabel[k]=Collabel[i];
         
         Cnonzeobestconfidence[k]=nonzeobestconfidence[i];
         foutlabel<<cRowlabel[k]<<" "<<cCollabel[k]<<endl;
         foutclass<<Cnonzerobestclass[k]<<endl;

         k=k+1;
cout<<"k:"<<k<<endl;
   //}
}


}

if (global==0)
{
///*
 for(int i=0;i<num;i++)
 {  
//if(Rowlabel[i]<winm/1.5&&Rowlabel[i]>winm*0.33&&Collabel[i]>winn/3&&Collabel[i]<winn*1) 


if (partloop==1)
{
cout<<"partloop:"<<partloop<<endl;
if (Rowlabel[i]<part1schp[3]&&Rowlabel[i]>part1schp[1]&&Collabel[i]>part1schp[0]&&Collabel[i]<part1schp[2])

       {
         Cnonzerobestclass[k]=nonzeobestclass[i];
//cout<<"Cnonzerobestclass"<<Cnonzerobestclass[k]<<endl;
         cRowlabel[k]=Rowlabel[i];
         cCollabel[k]=Collabel[i];
         
         Cnonzeobestconfidence[k]=nonzeobestconfidence[i];
         foutlabel<<cRowlabel[k]<<" "<<cCollabel[k]<<endl;
         foutclass<<Cnonzerobestclass[k]<<endl;

         k=k+1;
cout<<"k:"<<k<<endl;
   }

}

if (partloop==2)
{

if (Rowlabel[i]<part2schp[3]&&Rowlabel[i]>part2schp[1]&&Collabel[i]>part2schp[0]&&Collabel[i]<part2schp[2])

       {
         Cnonzerobestclass[k]=nonzeobestclass[i];
//cout<<"Cnonzerobestclass"<<Cnonzerobestclass[k]<<endl;
         cRowlabel[k]=Rowlabel[i];
         cCollabel[k]=Collabel[i];
         
         Cnonzeobestconfidence[k]=nonzeobestconfidence[i];
         foutlabel<<cRowlabel[k]<<" "<<cCollabel[k]<<endl;
         foutclass<<Cnonzerobestclass[k]<<endl;

         k=k+1;


   }

}

if (partloop==3)
{

if (Rowlabel[i]<part3schp[3]&&Rowlabel[i]>part3schp[1]&&Collabel[i]>part3schp[0]&&Collabel[i]<part3schp[2])

       {
         Cnonzerobestclass[k]=nonzeobestclass[i];
//cout<<"Cnonzerobestclass"<<Cnonzerobestclass[k]<<endl;
         cRowlabel[k]=Rowlabel[i];
         cCollabel[k]=Collabel[i];
         
         Cnonzeobestconfidence[k]=nonzeobestconfidence[i];
         foutlabel<<cRowlabel[k]<<" "<<cCollabel[k]<<endl;
         foutclass<<Cnonzerobestclass[k]<<endl;

         k=k+1;
   }

}


 }
//*/
}





matrix<float> Cnonzeroconfidence;
Cnonzeroconfidence.resize(k,n);
k=0;


if (global==1)
{

 for(int i=0;i<num;i++)
 {  
//if (Rowlabel[i]<winm/1.5&&Rowlabel[i]>winm*0.33&&Collabel[i]>winn/3&&Collabel[i]<winn*1) 

       //{
        for (int j=0;j<n;j++)
{
Cnonzeroconfidence(k,j)=nonzeroconfidence(i,j);
//cout<<"Cnonzeroconfidence:"<<Cnonzeroconfidence(k,j)<<endl;
//cout<<"n:"<<n<<endl;
}
         k=k+1;
   }
//}


}



if (global==0)
{
///*
 for(int i=0;i<num;i++)
 {  


if (partloop==1)
{

if (Rowlabel[i]<part1schp[3]&&Rowlabel[i]>part1schp[1]&&Collabel[i]>part1schp[0]&&Collabel[i]<part1schp[2])

       {
        for (int j=0;j<n;j++)
{
Cnonzeroconfidence(k,j)=nonzeroconfidence(i,j);
//cout<<"Cnonzeroconfidence:"<<Cnonzeroconfidence(k,j)<<endl;
//cout<<"n:"<<n<<endl;
}
         k=k+1;
   }
}


if (partloop==2)
{

if (Rowlabel[i]<part2schp[3]&&Rowlabel[i]>part2schp[1]&&Collabel[i]>part2schp[0]&&Collabel[i]<part2schp[2])

       {
        for (int j=0;j<n;j++)
{
Cnonzeroconfidence(k,j)=nonzeroconfidence(i,j);
//cout<<"Cnonzeroconfidence:"<<Cnonzeroconfidence(k,j)<<endl;
//cout<<"n:"<<n<<endl;
}
         k=k+1;
   }
}

if (partloop==3)
{

if (Rowlabel[i]<part3schp[3]&&Rowlabel[i]>part3schp[1]&&Collabel[i]>part3schp[0]&&Collabel[i]<part3schp[2])

       {
        for (int j=0;j<n;j++)
{
Cnonzeroconfidence(k,j)=nonzeroconfidence(i,j);
//cout<<"Cnonzeroconfidence:"<<Cnonzeroconfidence(k,j)<<endl;
//cout<<"n:"<<n<<endl;
}
         k=k+1;
   }
}


 }
//*/
}



//////compute aveconf1
 float evaconf1[200];
//float aveconf1[200];
float sum=0;
for(int j=0;j<n-1;j++)
 { evaconf1[j]=0;
 for(int i=0;i<k;i++)
 {  evaconf1[j]=evaconf1[j]+Cnonzeroconfidence(i,j+1);
}
}

for(int j=0;j<n-1;j++)
 { sum=sum+evaconf1[j];
}

for(int j=0;j<n-1;j++)
 { 

if (sum!=0)
{
evaconf1[j]=evaconf1[j]/sum;
}
else
{
evaconf1[j]=0;
}

cout<<"evaconf1:"<<evaconf1[j]<<endl;
}






foutlabel.close();
foutclass.close();
//////compute aveconf
 float preclassnum[200];
 //std::vector<float> aveconf;
 //std::vector<float> eva;
//std::vector<float> evaconf;
float aveconf[200];
 for(int i=0;i<200;i++)
{
      aveconf[i]=0;
      preclassnum[i]=0;
}
float eva[k];
extern float evaconf[200];
//cout<<"k"<<k<<endl;
for(int i=0;i<k;i++)
{
    preclassnum[Cnonzerobestclass[i]-1]=preclassnum[Cnonzerobestclass[i]-1]+1;

    aveconf[Cnonzerobestclass[i]-1]=aveconf[Cnonzerobestclass[i]-1]+Cnonzeobestconfidence[i];
//cout<<Cnonzeobestconfidence[i]<<endl<<aveconf[Cnonzerobestclass[i]]<<endl;

}


for(int i=0;i<n-1;i++)
{
    if(aveconf[i]!=0)
{//cout<<"aveconf"<<aveconf[i]<<endl;
aveconf[i]=aveconf[i]/preclassnum[i];
//cout<<"i"<<i<<endl<<preclassnum[i]<<endl;
    }
}

float evasum=0;
for(int i=0;i<n-1;i++)
{
eva[i]=aveconf[i]*preclassnum[i];

evasum=evasum+eva[i];

}

int maxnum;
float max;
if (evasum!=0)
{
for(int i=0;i<n-1;i++)
{
evaconf[i]=eva[i]/evasum;


if (i==0)
{max=evaconf[0];
maxnum=1;
}
else if(evaconf[i]>max)
{max=evaconf[i];
maxnum=i+1;
}
cout<<i+1<<" "<<evaconf[i]<<endl;


}
}
else
{
for(int i=0;i<n-1;i++)
{
evaconf[i]=0;

}
maxnum=0;
max=0;
}


cout<<"Car single part Result:"<<maxnum<<" "<<"confidence:"<<max<<endl;
/*
*/
//return evaconf[18];

}

