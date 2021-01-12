

#include "nodehyperplane.h"
#include "nodegini.h"
#include "nodeinfogain.h"
#include "utilities.h"
#include <boost/foreach.hpp>
#include <algorithm>
#ifdef WIN32
inline double round(double x)
{
    return (x-floor(x))>0.5 ? ceil(x) : floor(x);
}
#endif
using namespace std;

NodeHyperPlane::NodeHyperPlane(const HyperParameters &hp, int depth) : Node(hp, depth), m_bestThreshold( 0.0 )
{
//	visited=false;
}

NodeHyperPlane::NodeHyperPlane(const HyperParameters &hp, int depth, int reset) : Node(hp, depth, reset), m_bestThreshold( 0.0 )
{
}

NodeHyperPlane::NodeHyperPlane(const HyperParameters &hp, int depth, float fbestThreshold, const std::vector<int> & vbestFeatures, const std::vector<float> & vbestWeights)
		: Node(hp, depth),
		m_bestThreshold( fbestThreshold ),
		m_bestFeatures(vbestFeatures),
		m_bestWeights(vbestWeights)
{
//	visited=false;

}

NodeHyperPlane::NodeHyperPlane(const HyperParameters &hp, int depth, int reset, float fbestThreshold, const std::vector<int> & vbestFeatures, const std::vector<float> & vbestWeights)
		: Node(hp, depth, reset),
		m_bestThreshold( fbestThreshold ),
		m_bestFeatures(vbestFeatures),
		m_bestWeights(vbestWeights)
{
//	visited=false;

}
/////////////////////////yunxing///////////////////////
NodeHyperPlane::NodeHyperPlane(const HyperParameters &hp, int reset, const xmlNodePtr nodeNode) : Node(hp,0,reset)
{
    m_isLeaf = (readStringProp(nodeNode,"isLeaf") == "true") ? true : false;
//    visited=false;
    if (m_isLeaf)
    {//cout<<"test1"<<endl;
        m_nodeLabel = readIntProp( nodeNode, "label", 0 );
        xmlNodePtr cur = nodeNode->xmlChildrenNode;
        while ( cur != 0 )
        {
            if ( xmlStrcmp( cur->name, reinterpret_cast<const xmlChar*>( "confidence" ) ) == 0 )
            {
                m_nodeConf.push_back( static_cast<float>(readDoubleProp( cur, "conf", 0 )) );
            }
            cur = cur->next;
        }
    }
    else
    {
        xmlNodePtr cur = nodeNode->xmlChildrenNode;
        while ( cur != 0 ){

			if ( xmlStrcmp( cur->name, reinterpret_cast<const xmlChar*>( "data" ) ) == 0 ){
				
				xmlNodePtr curChild = cur->xmlChildrenNode;
				while ( curChild != 0 ){
					if ( xmlStrcmp( curChild->name, reinterpret_cast<const xmlChar*>( "bestfeatures" ) ) == 0 ){
						//cur = cur->next;
						xmlNodePtr curGrandChild = curChild->xmlChildrenNode;
						while ( curGrandChild != 0 ){
							if ( xmlStrcmp( curGrandChild->name, reinterpret_cast<const xmlChar*>( "feature" ) ) == 0 )
							{
								m_bestFeatures.push_back( static_cast<int>(readDoubleProp( curGrandChild, "value", 0 )) );
							}
							curGrandChild = curGrandChild->next;
						}
					}
					else if( xmlStrcmp( curChild->name, reinterpret_cast<const xmlChar*>( "bestweights" ) ) == 0 ){
						//cur = cur->next;
						
						xmlNodePtr curGrandChild = curChild->xmlChildrenNode;
						while ( curGrandChild != 0 ){
							if ( xmlStrcmp( curGrandChild->name, reinterpret_cast<const xmlChar*>( "weight" ) ) == 0 )
							{
								m_bestWeights.push_back( static_cast<float>(readDoubleProp( curGrandChild, "value", 0 )) );
							}
							curGrandChild = curGrandChild->next;
						}
					}
					else if( xmlStrcmp( curChild->name, reinterpret_cast<const xmlChar*>( "bestthreshold" ) ) == 0 ){
						float fbestThreshold = static_cast<float>(readDoubleProp( curChild, "value", 0 ));
						m_bestThreshold = fbestThreshold;
					}

					curChild = curChild->next;
				}

			}
            else if ( xmlStrcmp( cur->name, reinterpret_cast<const xmlChar*>( "node" ) ) == 0 )
            {
                const std::string childNode = readStringProp(cur,"child");
                if ( childNode == "left" )
                {
                    const std::string type = readStringProp(cur,"type");
                    if (type == NODE_GINI)
                    {
                        m_leftChildNode = NodeGini::Ptr(new NodeGini(m_hp,-1,cur));
                    }
                    else if (type == NODE_INFO_GAIN)
                    {
                        m_leftChildNode = NodeInfoGain::Ptr(new NodeInfoGain(m_hp,-1, cur));
                    }
                    else if (type == NODE_HYPER_PLANE)
                    {
                        m_leftChildNode = NodeHyperPlane::Ptr(new NodeHyperPlane(m_hp, -1, cur));
                    }
                }
                else
                {
                    const std::string type = readStringProp(nodeNode,"type");
                    if (type == NODE_GINI)
                    {
                        m_rightChildNode = NodeGini::Ptr(new NodeGini(m_hp,-1,cur));
                    }
                    else if (type == NODE_INFO_GAIN)
                    {
                        m_rightChildNode = NodeInfoGain::Ptr(new NodeInfoGain(m_hp,-1,cur));
                    }
                    else if (type == NODE_HYPER_PLANE)
                    {
                        m_rightChildNode = NodeHyperPlane::Ptr(new NodeHyperPlane(m_hp, -1, cur));
                    }
                }
            }
            cur = cur->next;
        }
    }
}

xmlNodePtr NodeHyperPlane::saveBestFeatures() const
{
    xmlNodePtr node = xmlNewNode( NULL, reinterpret_cast<const xmlChar*>( "bestfeatures" ) );

	std::vector<int>::const_iterator it(this->m_bestFeatures.begin()), end(this->m_bestFeatures.end());

	for(; it != end; it++){
		xmlNodePtr childNode = xmlNewNode( NULL, reinterpret_cast<const xmlChar*>( "feature" ) );
		addIntProp(childNode, "value", *it);
		xmlAddChild(node, childNode);
	}

    return node;
}

xmlNodePtr NodeHyperPlane::saveBestWeights() const {

    xmlNodePtr node = xmlNewNode( NULL, reinterpret_cast<const xmlChar*>( "bestweights" ) );

	std::vector<float>::const_iterator it(this->m_bestWeights.begin()), end(this->m_bestWeights.end());

	for(; it != end; it++){
		xmlNodePtr childNode = xmlNewNode( NULL, reinterpret_cast<const xmlChar*>( "weight" ) );
		addDoubleProp(childNode, "value", static_cast<float>(*it));
		xmlAddChild(node, childNode);
	}

    return node;
}

xmlNodePtr NodeHyperPlane::saveBestTreshold() const {

    xmlNodePtr node = xmlNewNode( NULL, reinterpret_cast<const xmlChar*>( "bestthreshold" ) );
	addDoubleProp(node, "value", static_cast<float>(this->m_bestThreshold));

    return node;
}

xmlNodePtr NodeHyperPlane::saveFeature() const {

	xmlNodePtr node = xmlNewNode( NULL, reinterpret_cast<const xmlChar*>( "data" ) );
	xmlAddChild(node, saveBestFeatures());
	xmlAddChild(node, saveBestWeights());
	xmlAddChild(node, saveBestTreshold());

    return node;
}


xmlNodePtr NodeHyperPlane::save() const
{
    xmlNodePtr node = xmlNewNode( NULL, reinterpret_cast<const xmlChar*>( "node" ) );
    xmlNewProp( node, reinterpret_cast<const xmlChar*>( "type" ), reinterpret_cast<const xmlChar*>( NODE_HYPER_PLANE ) );
    const char* isLeaf = (m_isLeaf) ? "true" : "false";
    xmlNewProp( node, reinterpret_cast<const xmlChar*>( "isLeaf" ),
                reinterpret_cast<const xmlChar*>( isLeaf ) );
    if (!m_isLeaf)
    {
        xmlAddChild(node, saveFeature());
        xmlNodePtr leftChildNode = m_leftChildNode->save();
        xmlNewProp( leftChildNode, reinterpret_cast<const xmlChar*>( "child" ),
                    reinterpret_cast<const xmlChar*>( LEFT_CHILD_NODE ) );
        xmlAddChild( node, leftChildNode );

        xmlNodePtr rightChildNode = m_rightChildNode->save();
        xmlNewProp( rightChildNode, reinterpret_cast<const xmlChar*>( "child" ),
                    reinterpret_cast<const xmlChar*>( RIGHT_CHILD_NODE ) );
        xmlAddChild( node, rightChildNode );
    }
    else
    {
        addIntProp( node, "label", m_nodeLabel);
        std::vector<float>::const_iterator it(m_nodeConf.begin()),end(m_nodeConf.end());
        int idx = 0;
        for (;it != end;it++,idx++)
        {
            xmlAddChild(node,saveConfidence(idx,*it));
        }
    }

    return node;
}


/*
xmlNodePtr NodeHyperPlane::save(FILE *fp) const
{
	typedef boost::shared_ptr<NodeHyperPlane> PtrHyperPlane;

	std::vector<PtrHyperPlane> stack_of_nodes;
	std::vector<bool> left_child;

	PtrHyperPlane child_node;
	bool child_bool;

	fprintf(fp,"\t\t<node type=\"nodeHyperPlane\" isLeaf=\"%s\">\n",(m_isLeaf) ? "true" : "false");

	fprintf(fp,"\t\t\t<data>\n");
	fprintf(fp,"\t\t\t\t<bestfeatures>\n");
	for (int i=0;i!=m_bestFeatures.size();i++)	{
		fprintf(fp,"\t\t\t\t\t<feature value=\"%d\"/>\n",m_bestFeatures[i]);
	}
	fprintf(fp,"\t\t\t\t</bestfeatures>\n");
	fprintf(fp,"\t\t\t\t<bestweights>\n");
	for (int i=0;i!=m_bestWeights.size();i++)	{
		fprintf(fp,"\t\t\t\t\t<weight value=\"%f\"/>\n",m_bestWeights[i]);
	}
	fprintf(fp,"\t\t\t\t</bestweights>\n");

	fprintf(fp,"\t\t\t\t<bestthreshold value=\"%f\"/>\n",m_bestThreshold);
	fprintf(fp,"\t\t\t</data>\n");

	child_node = boost::dynamic_pointer_cast<NodeHyperPlane>(m_leftChildNode);
	child_bool = true;

	while (child_node->m_isLeaf == false)	{
		stack_of_nodes.push_back(child_node);
		left_child.push_back(true);

		if (child_node->m_isLeaf == false) {
			fprintf(fp,"\t\t<node type=\"nodeHyperPlane\" isLeaf=\"%s\" child=\"%s\">\n",(child_node->m_isLeaf) ? "true" : "false",(child_bool) ? "left" : "right");

			fprintf(fp,"\t\t\t<data>\n");
			fprintf(fp,"\t\t\t\t<bestfeatures>\n");
			for (int i=0;i!=child_node->m_bestFeatures.size();i++)	{
				fprintf(fp,"\t\t\t\t\t<feature value=\"%d\"/>\n",child_node->m_bestFeatures[i]);
			}
			fprintf(fp,"\t\t\t\t</bestfeatures>\n");
			fprintf(fp,"\t\t\t\t<bestweights>\n");
			for (int i=0;i!=child_node->m_bestWeights.size();i++)	{
				fprintf(fp,"\t\t\t\t\t<weight value=\"%f\"/>\n",child_node->m_bestWeights[i]);
			}
			fprintf(fp,"\t\t\t\t</bestweights>\n");

			fprintf(fp,"\t\t\t\t<bestthreshold value=\"%f\"/>\n",child_node->m_bestThreshold);
			fprintf(fp,"\t\t\t</data>\n");
			
		}
		else	{
		}


		child_node = boost::dynamic_pointer_cast<NodeHyperPlane>(child_node->m_leftChildNode);
		child_bool = true;
	}

	stack_of_nodes.push_back(child_node);
	left_child.push_back(true);

	while (stack_of_nodes.size() > 0)	{
		child_node = stack_of_nodes.back();
		child_bool = left_child.back();

		stack_of_nodes.pop_back();
		left_child.pop_back();

		if (child_node->m_isLeaf == false) {
			fprintf(fp,"\t\t</node>\n");

			child_node = boost::dynamic_pointer_cast<NodeHyperPlane>(child_node->m_rightChildNode);
			child_bool = false;

			while (child_node->m_isLeaf == false)	{
				stack_of_nodes.push_back(child_node);
				left_child.push_back(child_bool);

				if (child_node->m_isLeaf == false) {
					fprintf(fp,"\t\t<node type=\"nodeHyperPlane\" isLeaf=\"%s\" child=\"%s\">\n",(child_node->m_isLeaf) ? "true" : "false",(child_bool) ? "left" : "right");

					fprintf(fp,"\t\t\t<data>\n");
					fprintf(fp,"\t\t\t\t<bestfeatures>\n");
					for (int i=0;i!=child_node->m_bestFeatures.size();i++)	{
						fprintf(fp,"\t\t\t\t\t<feature value=\"%d\"/>\n",child_node->m_bestFeatures[i]);
					}
					fprintf(fp,"\t\t\t\t</bestweights>\n");
					for (int i=0;i!=child_node->m_bestWeights.size();i++)	{
						fprintf(fp,"\t\t\t\t\t<weight value=\"%f\"/>\n",child_node->m_bestWeights[i]);
					}
					fprintf(fp,"\t\t\t\t</bestweights>\n");

					fprintf(fp,"\t\t\t\t<bestthreshold value=\"%f\"/>\n",child_node->m_bestThreshold);
					fprintf(fp,"\t\t\t</data>\n");
				}
				child_node = boost::dynamic_pointer_cast<NodeHyperPlane>(child_node->m_leftChildNode);
				child_bool = true;
			}
			stack_of_nodes.push_back(child_node);
			left_child.push_back(true);
		}
		else	{
			fprintf(fp,"\t\t<node type=\"nodeHyperPlane\" isLeaf=\"%s\" label=\"%d\" child=\"%s\">\n",(child_node->m_isLeaf) ? "true" : "false",child_node->m_nodeLabel,(child_bool) ? "left" : "right");
			for (int i=0;i!=child_node->m_nodeConf.size();i++)	{
				fprintf(fp,"\t\t\t<confidence class=\"%d\" conf=\"%f\"/>\n",i,child_node->m_nodeConf[i]);
			}
			fprintf(fp,"\t\t</node>\n");
		}

	}


	fprintf(fp,"\t\t</node>\n");
*/


xmlNodePtr NodeHyperPlane::save(FILE *fp) const
{
	// written by ZEESHAN
	// tree traversal is done in PREORDER
	// the </node> are printed by traversing the same tree in POSTORDER - when a leaf is reached in preorder, advance postorder algo (from wikipedia) until the last entry of postorder queue
	// is another leaf
	
	// visited flag based postorder implementation - POSTORDER
/*	iterativePostorder(rootNode)
	  nodeStack.push(rootNode)
	  while (! nodeStack.empty())
	    currNode = nodeStack.peek()
	    if ((currNode.left != null) and (currNode.left.visited == false))
	      nodeStack.push(currNode.left)
	    else 
	      if ((currNode.right != null) and (currNode.right.visited == false))
		nodeStack.push(currNode.right)
	      else
		print currNode.value
		currNode.visited := true
		nodeStack.pop()	*/

	typedef boost::shared_ptr<NodeHyperPlane> PtrHyperPlane;

	std::vector<PtrHyperPlane> stack_of_nodes;
	std::vector<bool> left_child;
//	std::vector<PtrHyperPlane> stack_postorder;
//	std::vector<bool> stack_whether_root;
//	std::vector<int> stack_nodeids;

	PtrHyperPlane child_node,lchild_node,rchild_node,rough_node;
	bool child_bool;

	int count_node_tag_starts,count_node_tag_ends;
	count_node_tag_starts = count_node_tag_ends = 0;

//	PtrHyperPlane current_postorder;//,lchild_postorder,rchild_postorder;
//	bool current_whether_root;

//	bool advance_postorder,advance_preorder;

	int id = 0;
	bool root=true;

	do 	{
		/////////////////////////
		// PREORDER
		/////////////////////////
//		if (advance_preorder == true)	{
		if (root==true)	{	// for the root node
			root=false;

			count_node_tag_starts++;
			fprintf(fp,"\t\t<node type=\"nodeHyperPlane\" isLeaf=\"%s\">\n",(m_isLeaf) ? "true" : "false");
			
			fprintf(fp,"\t\t\t<data>\n");
			fprintf(fp,"\t\t\t\t<bestfeatures>\n");
			for (int i=0;i!=m_bestFeatures.size();i++)	{
				fprintf(fp,"\t\t\t\t\t<feature value=\"%d\"/>\n",m_bestFeatures[i]);
			}
			fprintf(fp,"\t\t\t\t</bestfeatures>\n");
			fprintf(fp,"\t\t\t\t<bestweights>\n");
			for (int i=0;i!=m_bestWeights.size();i++)	{
				fprintf(fp,"\t\t\t\t\t<weight value=\"%.12e\"/>\n",m_bestWeights[i]);
			}
			fprintf(fp,"\t\t\t\t</bestweights>\n");

			fprintf(fp,"\t\t\t\t<bestthreshold value=\"%.12e\"/>\n",m_bestThreshold);
			fprintf(fp,"\t\t\t</data>\n");

			child_node = boost::dynamic_pointer_cast<NodeHyperPlane>(m_rightChildNode);
			child_bool = false;
			child_node->m_parent = boost::dynamic_pointer_cast<NodeHyperPlane>(child_node);
			child_node->m_child_of_root = true; 
			child_node->m_left_child = false;

			stack_of_nodes.push_back(child_node);
			left_child.push_back(child_bool);

			child_node = boost::dynamic_pointer_cast<NodeHyperPlane>(m_leftChildNode);
			child_bool = true;
			child_node->m_parent = boost::dynamic_pointer_cast<NodeHyperPlane>(child_node);
			child_node->m_child_of_root = true; 
			child_node->m_left_child = true;

			stack_of_nodes.push_back(child_node);
			left_child.push_back(child_bool);

		}
		else	{		// for all non-root nodes
			child_node = stack_of_nodes.back();
			child_bool = left_child.back();

			stack_of_nodes.pop_back();
			left_child.pop_back();

			if (child_node->m_isLeaf == false) {
				count_node_tag_starts++;
				fprintf(fp,"\t\t<node type=\"nodeHyperPlane\" isLeaf=\"%s\" child=\"%s\">\n",(child_node->m_isLeaf) ? "true" : "false",(child_bool) ? "left" : "right");

				fprintf(fp,"\t\t\t<data>\n");
				fprintf(fp,"\t\t\t\t<bestfeatures>\n");
				for (int i=0;i!=child_node->m_bestFeatures.size();i++)	{
					fprintf(fp,"\t\t\t\t\t<feature value=\"%d\"/>\n",child_node->m_bestFeatures[i]);
				}
				fprintf(fp,"\t\t\t\t</bestfeatures>\n");
				fprintf(fp,"\t\t\t\t<bestweights>\n");
				for (int i=0;i!=child_node->m_bestWeights.size();i++)	{
					fprintf(fp,"\t\t\t\t\t<weight value=\"%.12e\"/>\n",child_node->m_bestWeights[i]);
				}
				fprintf(fp,"\t\t\t\t</bestweights>\n");

				fprintf(fp,"\t\t\t\t<bestthreshold value=\"%.12e\"/>\n",child_node->m_bestThreshold);
				fprintf(fp,"\t\t\t</data>\n");

				rchild_node = boost::dynamic_pointer_cast<NodeHyperPlane>(child_node->m_rightChildNode);
				child_bool = false;
				rchild_node->m_parent = boost::dynamic_pointer_cast<NodeHyperPlane>(child_node);
				rchild_node->m_child_of_root = false; 
				rchild_node->m_left_child = false;

				stack_of_nodes.push_back(rchild_node);
				left_child.push_back(child_bool);

				lchild_node = boost::dynamic_pointer_cast<NodeHyperPlane>(child_node->m_leftChildNode);
				child_bool = true;
				lchild_node->m_parent = boost::dynamic_pointer_cast<NodeHyperPlane>(child_node); 
				lchild_node->m_child_of_root = false; 
				lchild_node->m_left_child = true;

				stack_of_nodes.push_back(lchild_node);
				left_child.push_back(child_bool);
			}
			else {
				count_node_tag_starts++;
				fprintf(fp,"\t\t<node type=\"nodeHyperPlane\" isLeaf=\"%s\" label=\"%d\" child=\"%s\">\n",(child_node->m_isLeaf) ? "true" : "false",child_node->m_nodeLabel,(child_bool) ? "left" : "right");
				for (int i=0;i!=child_node->m_nodeConf.size();i++)	{
					fprintf(fp,"\t\t\t<confidence class=\"%d\" conf=\"%.12e\"/>\n",i,child_node->m_nodeConf[i]);
				}
				fprintf(fp,"\t\t</node>\n");
				count_node_tag_ends++;
				
				if (child_node->m_left_child == false)	{
					rough_node = child_node;
					do {
						rough_node = boost::dynamic_pointer_cast<NodeHyperPlane>(rough_node->m_parent);
						fprintf(fp,"\t\t</node>\n");
						count_node_tag_ends++;
					} while (rough_node->m_left_child==false && rough_node->m_child_of_root!=true);

				}
			}
		}
	}
	while (stack_of_nodes.size() > 0);

	for (int i=0;i!=count_node_tag_starts-count_node_tag_ends;i++)	{
		fprintf(fp,"\t\t</node>\n");
	}


}


std::pair<float, float> NodeHyperPlane::calcGiniAndThreshold(const std::vector<int>& labels,
        const std::vector<std::pair<float, int> >& responses)
{
    // Initialize the counters: left takes all at the begining
    double DGini, LGini, RGini, LTotal, RTotal, bestW0 = 0, bestDGini = 1e10;
    std::vector<double> LCount(m_hp.numClasses, 0.0), RCount(m_hp.numClasses, 0.0);

    if (m_hp.isExtreme)
    {
        std::vector<std::pair<float, int> >::const_iterator resIt(responses.begin()), resEnd(responses.end());
        bestW0 = 2.0 * randomDouble( 1.0 ) - 1.0;
        RTotal = 0;
        LTotal = 0;
        for (; resIt != resEnd; resIt++)
        {
            if (resIt->first > bestW0)
            {
                RTotal++;
                RCount[labels[resIt->second]]++;
            }
            else
            {
                LTotal++;
                LCount[labels[resIt->second]]++;
            }
        }

        LGini = 0;
        RGini = 0;
        std::vector<double>::iterator LIt = LCount.begin(), RIt = RCount.begin(), end = LCount.end(), REnd = RCount.end();
        for (; LIt != end; LIt++, RIt++)      // Calculate Gini index
        {
            if (LTotal)
            {
                LGini += (*LIt/LTotal)*(1 - *LIt/LTotal);
            }
            if (RTotal)
            {
                RGini += (*RIt/RTotal)*(1 - *RIt/RTotal);
            }
        }

        bestDGini = (LTotal*LGini + RTotal*RGini)/responses.size();
    }
    else
    {
        RTotal = responses.size();
        LTotal = 0;

        // Count the number of samples in each class
        std::vector<std::pair<float, int> >::const_iterator resIt(responses.begin()), resEnd(responses.end()), tmpResIt;
        for (; resIt != resEnd; resIt++)
        {
            RCount[labels[resIt->second]]++;
        }

        // Loop over the sorted values and find the min DGini
        std::vector<double>::iterator LIt = LCount.begin(), RIt = RCount.begin(), end = LCount.end(), REnd = RCount.end();
        resIt = responses.begin();
        ++resIt;
        for (; resIt != resEnd; resIt++)
        {
            tmpResIt = resIt;
            --tmpResIt;

            RTotal--;
            LTotal++;
            RCount[labels[tmpResIt->second]]--;
            LCount[labels[tmpResIt->second]]++;

            if (resIt->first != tmpResIt->first)
            {
                LGini = 0;
                RGini = 0;
                LIt = LCount.begin();
                RIt = RCount.begin();
                for (; LIt != end; LIt++, RIt++)      // Calculate Gini index
                {
                    LGini += (*LIt/LTotal)*(1 - *LIt/LTotal);
                    RGini += (*RIt/RTotal)*(1 - *RIt/RTotal);
                }

                DGini = (LTotal*LGini + RTotal*RGini)/responses.size();
                if (DGini < bestDGini)
                {
                    bestDGini = DGini;
                    bestW0 = (resIt->first + tmpResIt->first)*0.5;
                }
            }
        }
    }

    return std::pair<float,float>((float)bestDGini,(float)bestW0);
}

std::pair<float, float> NodeHyperPlane::calcInfoGainAndThreshold(const std::vector<int>& labels,
        const std::vector<std::pair<float, int> >& responses)
{
    // Initialize the counters: left takes all at the begining
    double DInfo, LInfo, RInfo, LTotal, RTotal, bestW0 = 0.0, bestDInfo = 1e10;
    std::vector<double> LCount(m_hp.numClasses, 0.0), RCount(m_hp.numClasses, 0.0);

    RTotal = responses.size();
    LTotal = 0.0;
    // Count the number of samples in each class
    std::vector<std::pair<float, int> >::const_iterator resIt(responses.begin()), resEnd(responses.end()), tmpResIt;
    for (; resIt != resEnd; resIt++)
    {
        RCount[labels[resIt->second]]++;
    }

    // Loop over the sorted values and find the max DInfo
    std::vector<double>::iterator LIt = LCount.begin(), RIt = RCount.begin(), end = LCount.end(), REnd = RCount.end();
    resIt = responses.begin();
    ++resIt;
    for (; resIt != resEnd; resIt++)
    {
        tmpResIt = resIt;
        --tmpResIt;

        RTotal--;
        LTotal++;
        RCount[labels[tmpResIt->second]]--;
        LCount[labels[tmpResIt->second]]++;

        if (resIt->first != tmpResIt->first)
        {
            LInfo = 0.0;
            RInfo = 0.0;
            LIt = LCount.begin();
            RIt = RCount.begin();
            for (; LIt != end; LIt++, RIt++)      // Calculate Info index
            {
                if (*LIt)
                {
                    LInfo -= (*LIt/LTotal)*log(*LIt/LTotal);
                }
                if (*RIt)
                {
                    RInfo -= (*RIt/RTotal)*log(*RIt/RTotal);
                }
            }

            DInfo = (LTotal*LInfo + RTotal*RInfo)/responses.size();
            if (DInfo < bestDInfo)
            {
                bestDInfo = DInfo;
                bestW0 = (resIt->first + tmpResIt->first)*0.5;
            }
        }
    }

    return std::pair<float,float>((float)bestDInfo,(float)bestW0);
}

std::pair<float, float> NodeHyperPlane::calcInfoGainAndThreshold(const std::vector<int>& labels, const std::vector<double>& weights,
        const std::vector<std::pair<float, int> >& responses)
{
    // Initialize the counters: left takes all at the begining
    double DInfo, LInfo, RInfo, LTotal, RTotal, bestW0 = 0.0, bestDInfo = 1e10;
    std::vector<double> LCount(m_hp.numClasses, 0.0), RCount(m_hp.numClasses, 0.0);

    RTotal = 0.0;
    LTotal = 0.0;
    // Count the number of samples in each class
    std::vector<std::pair<float, int> >::const_iterator resIt(responses.begin()), resEnd(responses.end()), tmpResIt;
    for (; resIt != resEnd; resIt++)
    {
        RCount[labels[resIt->second]] += weights[resIt->second];
        RTotal += weights[resIt->second];
    }

    // Loop over the sorted values and find the max DInfo
    std::vector<double>::iterator LIt = LCount.begin(), RIt = RCount.begin(), end = LCount.end(), REnd = RCount.end();
    resIt = responses.begin();
    ++resIt;
    for (; resIt != resEnd; resIt++)
    {
        tmpResIt = resIt;
        --tmpResIt;

        RTotal -= weights[tmpResIt->second];
        LTotal += weights[tmpResIt->second];
        RCount[labels[tmpResIt->second]] -= weights[tmpResIt->second];
        LCount[labels[tmpResIt->second]] += weights[tmpResIt->second];

        if (resIt->first != tmpResIt->first)
        {
            LInfo = 0.0;
            RInfo = 0.0;
            LIt = LCount.begin();
            RIt = RCount.begin();
            for (; LIt != end; LIt++, RIt++)      // Calculate Info index
            {
                if (*LIt)
                {
                    LInfo -= (*LIt/LTotal)*log(*LIt/LTotal);
                }
                if (*RIt)
                {
                    RInfo -= (*RIt/RTotal)*log(*RIt/RTotal);
                }
            }

            DInfo = (LTotal*LInfo + RTotal*RInfo)/(LTotal + RTotal);
            if (DInfo < bestDInfo)
            {
                bestDInfo = DInfo;
                bestW0 = (resIt->first + tmpResIt->first)*0.5;
            }
        }
    }

    return std::pair<float,float>((float)bestDInfo,(float)bestW0);
}


std::pair<float, float> NodeHyperPlane::calcGiniAndThreshold(const std::vector<int>& labels, const std::vector<double>& weights,
        const std::vector<std::pair<float, int> >& responses, const bool useUnlabeledData)
{
    // Initialize the counters: left takes all at the begining
    double DGini, LGini, RGini, LTotal, RTotal, bestW0 = 0, bestDGini = 1e10;
    std::vector<double> LCount(m_hp.numClasses, 0.0), RCount(m_hp.numClasses, 0.0);

    RTotal = 0;
    LTotal = 0;

    // Count the number of samples in each class
    std::vector<std::pair<float, int> >::const_iterator resIt(responses.begin()), resEnd(responses.end()), tmpResIt;
    for (; resIt != resEnd; resIt++)
    {
        if (useUnlabeledData || resIt->second < m_hp.numLabeled)
        {
			/*
			CODIGO ORIGINAL COMENTADO.
            RCount[labels[resIt->second]] += weights[resIt->second];
            RTotal += weights[resIt->second];
			*/
			RCount[labels[resIt->second]] += weights[labels[resIt->second]];
            RTotal += weights[labels[resIt->second]];
        }
    }

    // Loop over the sorted values and find the min DGini
    std::vector<double>::iterator LIt = LCount.begin(), RIt = RCount.begin(), end = LCount.end(), REnd = RCount.end();
    resIt = responses.begin();
    ++resIt;
    for (; resIt != resEnd; resIt++)
    {
        if (useUnlabeledData || resIt->second < m_hp.numLabeled)
        {
            tmpResIt = resIt;
            --tmpResIt;

			/*
			CODIGO ORIGINAL COMENTADO
            RTotal -= weights[tmpResIt->second];
            LTotal += weights[tmpResIt->second];
            RCount[labels[tmpResIt->second]] -= weights[tmpResIt->second];
            LCount[labels[tmpResIt->second]] += weights[tmpResIt->second];
			*/
            RTotal -= weights[labels[tmpResIt->second]];
            LTotal += weights[labels[tmpResIt->second]];
            RCount[labels[tmpResIt->second]] -= weights[labels[tmpResIt->second]];
            LCount[labels[tmpResIt->second]] += weights[labels[tmpResIt->second]];

            if (resIt->first != tmpResIt->first)
            {
                LGini = 0;
                RGini = 0;
                LIt = LCount.begin();
                RIt = RCount.begin();
                for (; LIt != end; LIt++, RIt++)      // Calculate Gini index
                {
                    LGini += (*LIt/LTotal)*(1 - *LIt/LTotal);
                    RGini += (*RIt/RTotal)*(1 - *RIt/RTotal);
                }

                DGini = (LTotal*LGini + RTotal*RGini)/(LTotal + RTotal);
                if (DGini < bestDGini)
                {
                    bestDGini = DGini;
                    bestW0 = (resIt->first + tmpResIt->first)*0.5;
                }
            }
        }
    }

    return std::pair<float,float>((float)bestDGini,(float)bestW0);
}


void NodeHyperPlane::findHypotheses(const matrix<float>& data, const std::vector<int>& labels,
                                    const std::vector<int>& inBagSamples, const std::vector<int>& randFeatures, int numTries)
{
    std::vector<double> gini(m_hp.numRandomFeatures), thresholds(m_hp.numRandomFeatures);
    std::vector<int>::const_iterator it(randFeatures.begin());
    std::vector<int>::const_iterator end(randFeatures.end());

    std::vector<int>::const_iterator bagIt;
    std::vector<int>::const_iterator bagEnd(inBagSamples.end());

    double bestDGini = 1e10, bestThreshold = 0;
    std::pair<float,float> curGiniThresh;
    std::vector<std::pair<float, int> > responses;
    std::vector<float> bestWeights(randFeatures.size(),0.0);
    std::vector<float> tmpWeights(randFeatures.size(),0.0);
    float tmp = 0.0;
    for ( int i = 0; i < numTries; i++)
    {
        fillWithRandomNumbers(tmpWeights);
        responses.clear();
        responses.reserve(inBagSamples.size());
        bagIt = inBagSamples.begin();
        while ( bagIt != bagEnd )
        {
            tmp = 0.0;
            int counter = 0;
            BOOST_FOREACH(int feat, randFeatures)
            {
                tmp += data(*bagIt,feat)*tmpWeights[counter];
                counter++;
            }

            responses.push_back(std::pair<float, int>(tmp,*bagIt));
            ++bagIt;
        }

        sort(responses.begin(), responses.end());

        if (m_hp.useInfoGain)
        {
            curGiniThresh = calcInfoGainAndThreshold(labels, responses);
        }
        else
        {
            curGiniThresh = calcGiniAndThreshold(labels, responses);
        }

        if (curGiniThresh.first < bestDGini)
        {
            bestDGini     = curGiniThresh.first;
            bestThreshold = curGiniThresh.second;
            bestWeights   = tmpWeights;
        }

		gini[i]			= curGiniThresh.first;
		thresholds[i]	= curGiniThresh.second;
    }

    m_bestWeights = bestWeights;
    m_bestFeatures = randFeatures;
    m_bestThreshold = (float) bestThreshold;
}

void NodeHyperPlane::findHypothesesLU(const matrix<float>& data, const std::vector<int>& labels,
                                      const std::vector<int>& inBagSamples, const std::vector<int>& randFeatures, int numTries)
{
    std::vector<double> gini(m_hp.numRandomFeatures), thresholds(m_hp.numRandomFeatures);
    std::vector<int>::const_iterator it(randFeatures.begin());
    std::vector<int>::const_iterator end(randFeatures.end());

    std::vector<int>::const_iterator bagIt;
    std::vector<int>::const_iterator bagEnd(inBagSamples.end());

    double bestDGini = 1e10, bestThreshold = 0;
    std::pair<float,float> curGiniThresh;
    std::vector<std::pair<float, int> > responses;
    std::vector<float> bestWeights(randFeatures.size(),0.0);
    std::vector<float> tmpWeights(randFeatures.size(),0.0);
    float tmp = 0.0;
    for ( int i = 0; i < numTries; i++)
    {
        fillWithRandomNumbers(tmpWeights);
        responses.clear();
        responses.reserve(inBagSamples.size());
        bagIt = inBagSamples.begin();
        while ( bagIt != bagEnd )
        {
            tmp = 0.0;
            int counter = 0;
            BOOST_FOREACH(int feat, randFeatures)
            {
                tmp += data(*bagIt,feat)*tmpWeights[counter];
                counter++;
            }

            responses.push_back(std::pair<float, int>(tmp,*bagIt));
            ++bagIt;
        }

        sort(responses.begin(), responses.end());

        if (inBagSamples[0] < m_hp.numLabeled)
        {
            curGiniThresh = calcGiniAndThreshold(labels, responses);
        }
        else
        {
            curGiniThresh = calcClusterScoreAndThreshold(data, inBagSamples, responses);
        }
        if (curGiniThresh.first < bestDGini)
        {
            bestDGini = curGiniThresh.first;
            bestThreshold = curGiniThresh.second;
            bestWeights = tmpWeights;
        }
    }

    m_bestWeights = bestWeights;
    m_bestFeatures = randFeatures;
    m_bestThreshold = (float) bestThreshold;
}

std::pair<float, float> NodeHyperPlane::calcClusterScoreAndThreshold(const matrix<float>& data, const std::vector<int>& inBagSamples,
        const std::vector<double>& weights,
        const std::vector<std::pair<float, int> >& responses)
{
    // Find the mid-point using responsess
    int numResponse = responses.size();
    int midPointIndex = (int) round(numResponse/2 - 1);
    float threshold = (responses[midPointIndex].first + responses[midPointIndex + 1].first)/2;

    // Calculate the weighted cluster center for left and right splits
    double LWeight = 0, RWeight = 0;
    std::vector<float> LCenter(data.size2(), 0.0), RCenter(data.size2(), 0.0);
    for (int m = 0; m < (int) data.size2(); m++)
    {
        for (int n = 0; n < midPointIndex; n++)
        {
            LCenter[m] += (float)weights[responses[n].second]*data(responses[n].second, m);
            if (m == 0)
            {
                LWeight += weights[responses[n].second];
            }
        }
        LCenter[m] /= (float)LWeight;

        for (int n = midPointIndex; n < (int) responses.size(); n++)
        {
            RCenter[m] += (float)weights[responses[n].second]*data(responses[n].second, m);
            if (m == 0)
            {
                RWeight += weights[responses[n].second];
            }
        }
        RCenter[m] /= (float)RWeight;
    }

    // Calculate the weighted distance from each point to the centers
    float LScore = 0, RScore = 0;
    for (int m = 0; m < (int) data.size2(); m++)
    {
        for (int n = 0; n < midPointIndex; n++)
        {
            LScore += weights[responses[n].second]*pow((double) (data(responses[n].second, m) - LCenter[m]), 2.0);
        }
        LScore /= (float)LWeight;
        for (int n = midPointIndex; n < (int) responses.size(); n++)
        {
            RScore += (float)weights[responses[n].second]*pow((double) (data(responses[n].second, m) - RCenter[m]), 2.0);
        }
        RScore /= (float)RWeight;
    }

    return std::pair<float,float>(0.5f*(LScore + RScore), threshold);
}

std::pair<float, float> NodeHyperPlane::calcClusterScoreAndThreshold(const matrix<float>& data, const std::vector<int>& inBagSamples,
        const std::vector<std::pair<float, int> >& responses)
{
    // Find the mid-point using responsess
    int numResponse = responses.size();
    int midPointIndex = (int) round(numResponse/2 - 1);
    float threshold = (responses[midPointIndex].first + responses[midPointIndex + 1].first)/2;

    // Calculate the weighted cluster center for left and right splits
    double LWeight = 0, RWeight = 0;
    std::vector<float> LCenter(data.size2(), 0.0), RCenter(data.size2(), 0.0);
    for (int m = 0; m < (int) data.size2(); m++)
    {
        for (int n = 0; n < midPointIndex; n++)
        {
            LCenter[m] += data(responses[n].second, m);
            if (m == 0)
            {
                LWeight++;
            }
        }
        LCenter[m] /= LWeight;

        for (int n = midPointIndex; n < (int) responses.size(); n++)
        {
            RCenter[m] += data(responses[n].second, m);
            if (m == 0)
            {
                RWeight++;
            }
        }
        RCenter[m] /= RWeight;
    }

    // Calculate the weighted distance from each point to the centers
    float LScore = 0, RScore = 0;
    for (int m = 0; m < (int) data.size2(); m++)
    {
        for (int n = 0; n < midPointIndex; n++)
        {
            LScore += pow((double) (data(responses[n].second, m) - LCenter[m]), 2.0);
        }
        LScore /= LWeight;
        for (int n = midPointIndex; n < (int) responses.size(); n++)
        {
            RScore += pow((double) (data(responses[n].second, m) - RCenter[m]), 2.0);
        }
        RScore /= RWeight;
    }

    return std::pair<float,float>(0.5*(LScore + RScore), threshold);
}

void NodeHyperPlane::findHypotheses(const matrix<float>& data, const std::vector<int>& labels,
                                    const std::vector<double>& weights,
                                    const std::vector<int>& inBagSamples, const std::vector<int>& randFeatures, int numTries)
{
    std::vector<double> gini(m_hp.numRandomFeatures), thresholds(m_hp.numRandomFeatures);
    std::vector<int>::const_iterator it(randFeatures.begin());
    std::vector<int>::const_iterator end(randFeatures.end());
    std::vector<int>::const_iterator bagIt;
    std::vector<int>::const_iterator bagEnd(inBagSamples.end());

    double bestDGini = 1e10, bestThreshold = 0.0;
    std::pair<float,float> curGiniThresh;
    std::vector<std::pair<float, int> > responses;

    std::vector<float> bestWeights(randFeatures.size(),0.0);
    std::vector<float> tmpWeights(randFeatures.size(),0.0);
    float tmp = 0.0;
    bool doClustering = clusterOrGini(), useUnlabeledData = true;
    for ( int i = 0; i < numTries; i++)
    {
        fillWithRandomNumbers(tmpWeights);
        responses.clear();
        responses.reserve(inBagSamples.size());
        bagIt = inBagSamples.begin();
        while ( bagIt != bagEnd )
        {
            tmp = 0.0;
            int counter = 0;
            BOOST_FOREACH(int feat, randFeatures)
            {
                tmp += data(*bagIt,feat)*tmpWeights[counter];
                counter++;
            }

            responses.push_back(std::pair<float, int>(tmp,*bagIt));
            ++bagIt;
        }
        sort(responses.begin(), responses.end());

        if (!doClustering)
        {
            if (m_hp.useInfoGain)
            {
                curGiniThresh = calcInfoGainAndThreshold(labels, weights, responses);
            }
            else
            {
                curGiniThresh = calcGiniAndThreshold(labels, weights, responses, useUnlabeledData);
            }
        }
        else
        {
            curGiniThresh = calcClusterScoreAndThreshold(data, inBagSamples, weights, responses);
        }

        if (curGiniThresh.first < bestDGini)
        {
            bestDGini = curGiniThresh.first;
            bestThreshold = curGiniThresh.second;
            bestWeights = tmpWeights;
        }
    }

    m_bestFeatures  = randFeatures;
    m_bestWeights   = bestWeights;
    m_bestThreshold = (float) bestThreshold;
}

NODE_TRAIN_STATUS NodeHyperPlane::trainLU(const matrix<float>& data, const std::vector<int>& labels,
        std::vector<int>& inBagSamples, matrix<float>& confidences, std::vector<int>& predictions)
{
    bool doSplit = shouldISplitLU(labels,inBagSamples);
    NODE_TRAIN_STATUS myTrainingStatus = IS_NOT_LEAF;

    if ( doSplit )
    {
        m_isLeaf = false;

        //train here the node: Select random features and evaluate them
        std::vector<int> randFeatures = randPerm(data.size2(), m_hp.numProjFeatures );
        int numTries = m_hp.numRandomFeatures;// * (m_depth+1);
        findHypothesesLU(data, labels, inBagSamples, randFeatures, numTries);
        if (m_hp.verbose)
        {
            cout << "Node #: " << m_nodeIndex;
            cout << " and the threshold is: " << m_bestThreshold << " at depth " << m_depth << endl;
        }

        // split the data
        std::vector<int> leftNodeSamples, rightNodeSamples;
        evalNode(data,inBagSamples,leftNodeSamples,rightNodeSamples);

        // pass them to the left and right child, respectively
        m_leftChildNode = Ptr(new NodeHyperPlane(m_hp,m_depth + 1));
        m_rightChildNode = Ptr(new NodeHyperPlane(m_hp,m_depth + 1));

        m_leftChildNode->train(data,labels,leftNodeSamples,confidences,predictions);
        m_rightChildNode->train(data,labels,rightNodeSamples,confidences,predictions);
    }
    else
    {
        if (m_hp.verbose)
        {
            cout << "Node #: " << m_nodeIndex << " is terminal, at depth " << m_depth << endl;
        }

        // calc confidence, labels, etc
        m_isLeaf = true;
        myTrainingStatus = IS_LEAF;
		m_nodeConf.clear();
        m_nodeConf.resize(m_hp.numClasses, 0.0);
        int numNodeLabeled = 0;
        BOOST_FOREACH(int n, inBagSamples)
        {
            if (n < m_hp.numLabeled)
            {
                m_nodeConf[labels[n]]++;
                numNodeLabeled++;
            }
        }

        int bestClass = 0, tmpN = 0;
        float bestConf = 0;
        std::vector<float>::iterator confItr = m_nodeConf.begin(), confEnd = m_nodeConf.end();
        for (; confItr != confEnd; confItr++)
        {
            if (numNodeLabeled)
            {
                *confItr /= numNodeLabeled;
                if (*confItr > bestConf)
                {
                    bestConf = *confItr;
                    bestClass = tmpN;
                }
                tmpN++;
            }
            else
            {
                *confItr = 1.0/m_hp.numClasses;
            }
        }
        m_nodeLabel = bestClass;
//cout<<"test2"<<endl;
        BOOST_FOREACH(int n, inBagSamples)
        {
            predictions[n] = m_nodeLabel;
            tmpN = 0;
            BOOST_FOREACH(float conf, m_nodeConf)
            {
                confidences(n, tmpN) = conf;
                tmpN++;
            }
        }
    }

    return myTrainingStatus;
}

bool NodeHyperPlane::clusterOrGini()
{
    return false;
}
/////////////////////////////////step6evalnode/////////////////////////////////


void NodeHyperPlane::evalNode(const matrix<float>& data, const std::vector<int>& inBagSamples,
                              std::vector<int>& leftNodeSamples, std::vector<int>& rightNodeSamples)
{//cout<<"test"<<endl;
    float tmp;
    BOOST_FOREACH(int n, inBagSamples)
    {
        tmp = 0.0;
        int counter = 0;
        BOOST_FOREACH(int feat, m_bestFeatures)
        {
            tmp += data(n,feat)*m_bestWeights[counter];
            counter++;
        }

        if (tmp > m_bestThreshold)
        {
            rightNodeSamples.push_back(n);
        }
        else
        {
            leftNodeSamples.push_back(n);
        }
    }
}

NODE_TRAIN_STATUS NodeHyperPlane::train(const matrix<float>& data, const std::vector<int>& labels,
                                        std::vector<int>& inBagSamples, matrix<float>& confidences, std::vector<int>& predictions)
{//cout<<"test"<<endl;
    bool doSplit = shouldISplit(labels,inBagSamples);
    NODE_TRAIN_STATUS myTrainingStatus = IS_NOT_LEAF;

    if ( doSplit )
    {//cout<<"test2"<<endl;
        m_isLeaf = false;

        //train here the node: Select random features and evaluate them
        std::vector<int> randFeatures = randPerm(data.size2(), m_hp.numProjFeatures );
        int numTries = m_hp.numRandomFeatures;// * (m_depth+1);
        findHypotheses(data, labels, inBagSamples, randFeatures, numTries);
        if (m_hp.verbose)
        {
            cout << "Node #: " << m_nodeIndex;
            cout << " and the threshold is: " << m_bestThreshold << " at depth " << m_depth << endl;
        }

        // split the data
        std::vector<int> leftNodeSamples, rightNodeSamples;
        evalNode(data,inBagSamples,leftNodeSamples,rightNodeSamples);

        // pass them to the left and right child, respectively
        m_leftChildNode = Ptr(new NodeHyperPlane(m_hp,m_depth + 1));
        m_rightChildNode = Ptr(new NodeHyperPlane(m_hp,m_depth + 1));

        NODE_TRAIN_STATUS leftChildStatus = m_leftChildNode->train(data,labels,leftNodeSamples,confidences,predictions);
        NODE_TRAIN_STATUS rightChildStatus= m_rightChildNode->train(data,labels,rightNodeSamples,confidences,predictions);

    }
    else
    {
        if (m_hp.verbose)
        {
            cout << "Node #: " << m_nodeIndex << " is terminal, at depth " << m_depth << endl;
        }

        // calc confidence, labels, etc
        m_isLeaf = true;
        myTrainingStatus = IS_LEAF;

		m_nodeConf.clear();
        m_nodeConf.resize(m_hp.numClasses, 0.0);

        BOOST_FOREACH(int n, inBagSamples)
        {
            m_nodeConf[labels[n]]++;
        }

        int bestClass = 0, tmpN = 0;
        float bestConf = 0;
        std::vector<float>::iterator confItr = m_nodeConf.begin(), confEnd = m_nodeConf.end();
        for (; confItr != confEnd; confItr++)
        {
            *confItr /= inBagSamples.size();
            if (*confItr > bestConf)
            {
                bestConf = *confItr;
                bestClass = tmpN;
            }
            tmpN++;
        }
        m_nodeLabel = bestClass;

        BOOST_FOREACH(int n, inBagSamples)
        {
            predictions[n] = m_nodeLabel;
            tmpN = 0;
            BOOST_FOREACH(float conf, m_nodeConf)
            {
                confidences(n, tmpN) = conf;
                tmpN++;
            }
        }
    }

    return myTrainingStatus;
}

NODE_TRAIN_STATUS NodeHyperPlane::train(const matrix<float>& data, const std::vector<int>& labels, const std::vector<double>& weights,
                                        std::vector<int>& inBagSamples, matrix<float>& confidences, std::vector<int>& predictions)
{//cout<<"test"<<endl;
    bool doSplit = shouldISplit(labels,inBagSamples);
    NODE_TRAIN_STATUS myTrainingStatus = IS_NOT_LEAF;

    if ( doSplit )
    {
        m_isLeaf = false;

        //train here the node: Select random features and evaluate them
        std::vector<int> randFeatures = randPerm(data.size2(),m_hp.numProjFeatures );
        int numTries = m_hp.numRandomFeatures * (m_depth+1);
        findHypotheses(data, labels, weights, inBagSamples, randFeatures,numTries);
        if (m_hp.verbose)
        {
            cout << "Node #: " << m_nodeIndex;
            cout << " and the threshold is: " << m_bestThreshold << " at depth " << m_depth << endl;
        }

        // split the data
        std::vector<int> leftNodeSamples, rightNodeSamples;
        evalNode(data,inBagSamples,leftNodeSamples,rightNodeSamples);

        // pass them to the left and right child, respectively
        m_leftChildNode = Ptr(new NodeHyperPlane(m_hp,m_depth + 1));
        m_rightChildNode = Ptr(new NodeHyperPlane(m_hp,m_depth + 1));

        NODE_TRAIN_STATUS leftChildStatus = m_leftChildNode->train(data,labels,weights,leftNodeSamples,confidences,predictions);
        NODE_TRAIN_STATUS rightChildStatus= m_rightChildNode->train(data,labels,weights,rightNodeSamples,confidences,predictions);

    }
    else
    {
        if (m_hp.verbose)
        {
            cout << "Node #: " << m_nodeIndex << " is terminal, at depth " << m_depth << endl;
        }

        // calc confidence, labels, etc
        m_isLeaf = true;
        myTrainingStatus = IS_LEAF;

		m_nodeConf.clear();
        m_nodeConf.resize(m_hp.numClasses, 0.0);

        double totalW = 0;
        BOOST_FOREACH(int n, inBagSamples)
        {
			/*
			CODIGO ORIGINAL COMENTADO
            m_nodeConf[labels[n]] += weights[n];
            totalW += weights[n];
			*/
            m_nodeConf[labels[n]] += weights[labels[n]];
            totalW += weights[labels[n]];
        }

        int bestClass = 0, tmpN = 0;
        float bestConf = 0;
        std::vector<float>::iterator confItr = m_nodeConf.begin(), confEnd = m_nodeConf.end();
        for (; confItr != confEnd; confItr++)
        {
            *confItr /= (totalW + 1e-10);
            if (*confItr > bestConf)
            {
                bestConf = *confItr;
                bestClass = tmpN;
            }
            tmpN++;
        }
        m_nodeLabel = bestClass;

        BOOST_FOREACH(int n, inBagSamples)
        {
            predictions[n] = m_nodeLabel;
            tmpN = 0;
            BOOST_FOREACH(float conf, m_nodeConf)
            {
                confidences(n, tmpN) = conf;
                tmpN++;
            }
        }
    }

    return myTrainingStatus;
}

/////////////////////////////////step5eval/////////////////////////////////

void NodeHyperPlane::eval(const matrix<float>& data, const std::vector<int>& sampleIndeces,
                          matrix<float>& confidences, std::vector<int>& predictions)
{//cout<<"1"<<endl;
//cout<<"m_isLeaf"<<m_isLeaf<<endl;
    if (m_isLeaf)
    {//cout<<"leaf"<<endl;
        // Make predictions and confidences
        int tmpN;
       
        BOOST_FOREACH( int n, sampleIndeces)
        {//cout<<"n:"<<n<<endl;
            predictions[n] = m_nodeLabel;//////////////////////
            tmpN = 0;
 //cout<<"m_nodeConf:"<<m_nodeConf[20]<<endl;
            BOOST_FOREACH(float conf, m_nodeConf)////////keypoint
            {
                confidences(n, tmpN) = conf;////////////////////
                tmpN++;
///cout<<"tmpN:"<<tmpN<<endl;//////////////////////only to 20
//cout<<"conf:"<<conf<<endl;
//cout<<"m_nodeLabel"<<m_nodeLabel<<endl;
            }
        }
    }
    else
    {//cout<<"test"<<endl;
        // split the data
        std::vector<int> leftNodeSamples, rightNodeSamples;
        evalNode(data,sampleIndeces,leftNodeSamples,rightNodeSamples);
//////////////compute left right nodesamples///Ptr m_leftChildNode
        m_leftChildNode->eval(data,leftNodeSamples,confidences,predictions);//?
        m_rightChildNode->eval(data,rightNodeSamples,confidences,predictions);
    }
}

void NodeHyperPlane::getPath(const matrix<float>& data, const std::vector<int>& sampleIndeces, std::vector<std::vector<int> >& path)
{
    BOOST_FOREACH(int n, sampleIndeces)
    {
        path[n].push_back(m_nodeIndex);
    }

    if (!m_isLeaf)
    {
        // split the data
        std::vector<int> leftNodeSamples, rightNodeSamples;
        evalNode(data,sampleIndeces,leftNodeSamples,rightNodeSamples);

        m_leftChildNode->getPath(data,leftNodeSamples,path);
        m_rightChildNode->getPath(data,rightNodeSamples,path);
    }
}

void NodeHyperPlane::refine(const matrix<float>& data, const std::vector<int>& labels,
                            std::vector<int>& samples, matrix<float>& confidences, std::vector<int>& predictions)
{
    if ( !m_isLeaf )
    {
        // split the data
        std::vector<int> leftNodeSamples, rightNodeSamples;
        evalNode(data,samples,leftNodeSamples,rightNodeSamples);

        m_leftChildNode->refine(data,labels,leftNodeSamples,confidences,predictions);
        m_rightChildNode->refine(data,labels,rightNodeSamples,confidences,predictions);
    }
    else
    {//cout<<"test"<<endl;
       // calc confidence, labels, etc
		m_nodeConf.clear();
        m_nodeConf.resize(m_hp.numClasses,0.0);
        BOOST_FOREACH(int n, samples)
        {//cout<<"test"<<endl;
            m_nodeConf[labels[n]]++;
        }

        int bestClass = 0, tmpN = 0;
        float bestConf = 0;
        std::vector<float>::iterator confItr = m_nodeConf.begin(), confEnd = m_nodeConf.end();
        for (; confItr != confEnd; confItr++)
        {
            *confItr /= samples.size();
            if (*confItr > bestConf)
            {
                bestConf = *confItr;
                bestClass = tmpN;
            }
            tmpN++;
        }
        m_nodeLabel = bestClass;
//cout<<"test"<<endl;
        BOOST_FOREACH(int n, samples)
        {
            predictions[n] = m_nodeLabel;
            tmpN = 0;
            BOOST_FOREACH(float conf, m_nodeConf)
            {
                confidences(n, tmpN) = conf;
                tmpN++;
            }
        }
    }
}

void NodeHyperPlane::refine(const matrix<float>& data, const std::vector<int>& labels, const std::vector<double>& weights,
                            std::vector<int>& samples, matrix<float>& confidences, std::vector<int>& predictions)
{//cout<<"test"<<endl;
    // calc confidence, labels, etc
    if (!m_isLeaf)
    {
        // split the data
        std::vector<int> leftNodeSamples, rightNodeSamples;
        evalNode(data,samples,leftNodeSamples,rightNodeSamples);

        m_leftChildNode->refine(data,labels,weights,leftNodeSamples,confidences,predictions);
        m_rightChildNode->refine(data,labels,weights,rightNodeSamples,confidences,predictions);
    }
    else
    {//cout<<"test"<<endl;
		m_nodeConf.clear();
        m_nodeConf.resize(m_hp.numClasses,0.0);
        double totalW = 0;
        BOOST_FOREACH(int n, samples)
        {
			/*
			CODIGO ORIGINAL COMENTADO.
            m_nodeConf[labels[n]] += weights[n];
            totalW += weights[n];
			*/
            m_nodeConf[labels[n]] += weights[labels[n]];
            totalW += weights[labels[n]];
        }

        int bestClass = 0, tmpN = 0;
        float bestConf = 0;
        std::vector<float>::iterator confItr = m_nodeConf.begin(), confEnd = m_nodeConf.end();
        for (; confItr != confEnd; confItr++)
        {
            *confItr /= (totalW + 1e-10);
            if (*confItr > bestConf)
            {
                bestConf = *confItr;
                bestClass = tmpN;
            }
            tmpN++;
        }
        m_nodeLabel = bestClass;

        BOOST_FOREACH(int n, samples)
        {
            predictions[n] = m_nodeLabel;
            tmpN = 0;
            BOOST_FOREACH(float conf, m_nodeConf)
            {
                confidences(n, tmpN) = conf;
                tmpN++;
            }
        }
    }
}
