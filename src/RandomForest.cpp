#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include<math.h>
#include "forest.h"
#include "pairforest.h"
#include "tree.h"
#include "data.h"
#include "utilities.h"
#include "hyperparameters.h"
#include <libconfig.h++>
#include <cstdlib>
#include <libxml/tree.h>
#include <libxml/parser.h>
#include<boost/foreach.hpp>

using namespace std;
using namespace libconfig;

int hogfun1(int m, int n, float *pp, float *feat)
{
	//cout << "test1" << endl;
	 float Ix[200][200], Ied[200][200], Iphase[200][200];
	 float Iy[200][200];// Ied[500][500], Iphase[500][500];
	 float img[200][200];
	 float Cell[9][1000][1000];
	 float tmpx[26][26], tmped[26][26], tmpphase[26][26], f[9], sum1;
	 float imginput[110][110];
	   sum1 = 0;
	//cout << "test2" << endl;

	for (int i = 0; i < m; i++)
	{

		for (int j = 0; j < n; j++)
		{
			img[i][j] = *(pp + 200 * i + j);




		}
	}

	//cout << *pp << endl;

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			img[i][j] = sqrt(img[i][j]);      //伽马校正
		}
	}

	//
	//下面是求边缘


	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			int x1, x2, y1, y2;
			x1 = 1;
			x2 = -1;

			y1 = 1;
			y2 = -1;
			if (i == 0)
			{
				x1 = 1;
				x2 = 0;

			}
			if (i == m - 1)
			{
				x1 = 0;
				x2 = -1;

			}
			if (j == 0)
			{
				y1 = 1;
				y2 = 0;
			}
			if (j == n - 1)
			{
				y1 = 0;
				y2 = -1;
				//cout << Ix << endl << Iy;
			}

			Ix[i][j] = img[i + x1][j] - img[i + x2][j];
			Iy[i][j] = img[i][j + y1] - img[i][j + y2];
		}
	}
	//	cout << Ix[0][8] << endl;
	//
	//
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{

			Ied[i][j] = sqrt(pow(Ix[i][j], 2) + pow(Iy[i][j], 2));              //边缘强度
			//Iphase[i][j] = Iy[i][j] / Ix[i][j];              //边缘斜率，有些为inf,-inf,nan，其中nan需要再处理一下
		}
	}
	//cout << Ied[0][8] << endl;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (Ix[i][j] == 0 && Iy[i][j] == 0)
				Iphase[i][j] = 0;            //边缘斜率，有些为inf, -inf, nan，其中nan需要再处理一下
			else
				Iphase[i][j] = Iy[i][j] / Ix[i][j];
		}
	}
	//cout << Iphase[0][7] << endl;
	////下面是求cell
	int step = 16;                //step*step个像素作为一个单元
	float orient = 9;               //方向直方图的方向个数
	float jiao = 360 / orient, ang;        //每个方向包含的角度数
	int ii = 0;
	int jj = 0;
	float sum = 0;



	for (int i = 1; i <= m - step; i = i + step)
	{
		//如果处理的m/step不是整数，最好是i=1:step:m-step
		ii = 0;
		for (int j = 1; j <= n - step; j = j + step)      //注释同上
		{
			sum = 0;
			for (int s = 0; s < step; s++)
			{
				for (int t = 0; t < step; t++)
				{
					tmpx[s][t] = Ix[i - 1 + s][j - 1 + t];
					tmped[s][t] = Ied[i - 1 + s][j - 1 + t];
					tmpphase[s][t] = Iphase[i - 1 + s][j - 1 + t];
					sum = sum + tmped[s][t];
				}
			}
			for (int s = 0; s < step; s++)
			{
				for (int t = 0; t < step; t++)
				{
					tmped[s][t] = tmped[s][t] / sum;        //局部边缘强度归一化
				}
			}

			float Hist[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
			for (int p = 0; p < step; p++)
			{
				for (int q = 0; q < step; q++)
				{


					ang = atan(tmpphase[p][q]);    //atan求的是[-90 90]度之间

					ang = fmod(ang * 180 / 3.1415926 + 360, 360);    //全部变正，-90变270

					if (tmpx[p][q]<0)             //根据x方向确定真正的角度
					{
						if (ang<90)              //如果是第一象限
							ang = ang + 180;        //移到第三象限

						if (ang>270)              //如果是第四象限
							ang = ang - 180;        //移到第二象限

					}
					ang = ang + 0.0000001;          //防止ang为0
					int tmp = ceil(ang / jiao);

					Hist[tmp - 1] = Hist[tmp - 1] + tmped[p][q];   //ceil向上取整，使用边缘强度加权

				}
			}
			//if (i == 0 && j == 0 )//&& p == 0 && q == 10)

			//cout << Hist[0] << endl;
			//	//Hist=Hist/sum(Hist);    //方向直方图归一化，这一步可以没有，因为是组成block以后再进行归一化就可以
			for (int p = 0; p < orient; p++)
			{

				Cell[p][ii][jj] = Hist[p];       //放入Cell中




			}
			//	//Cell{ii,jj}=Hist;
			ii = ii + 1;                //针对Cell的y坐标循环变量
			//if (ii == 1 && jj == 1)//&& p == 0 && q == 10)

			//	cout <<"cell"<< Cell[8][0][1] << endl;
		}
		//针对Cell的x坐标循环变量
		jj = jj + 1;
	}


	//下面是求feature,2*2个cell合成一个block,没有显式的求block


	for (int i = 0; i < ii; i++)
	{
		for (int j = 0; j < jj; j++)
		{
			sum1 = 0;
			for (int t = 0; t < 1; t++)
			{
				for (int s = 0; s < 9; s++)
				{
					if (t == 0)
						f[t * 9 + s] = Cell[s][i][j];
					if (t == 1)
						f[t * 9 + s] = Cell[s][i][j + 1];
					if (t == 2)
						f[t * 9 + s] = Cell[s][i + 1][j];
					if (t == 3)
					{
						f[t * 9 + s] = Cell[s][i + 1][j + 1];


					}
					sum1 = sum1 + f[t * 9 + s];
				}
			}

			for (int i = 0; i < 9; i++)
			{
				f[i] = f[i] / sum1;//归一化
			}

			for (int t = 0; t < 9; t++)
			{

				//	feature[(i *(jj - 1) + j) * 36 + t] = f[t];
				*(feat + (i *(jj) + j) * 9 + t) = f[t];
			}



		}
	}



	//ofstream fout;
	//fout.open("E://feature.txt");

	//for (int in = 0; in < (ii - 1)*(jj - 1) * 36; in++)
	//{

	//fout << feature[in];

	//fout << endl;
	//if (in == 0)
	//cout << feature[in]<<endl;
	//}
	//fout.close();
	//int num = (ii - 1)*(jj - 1) * 9;
int num = (ii)*(jj) * 9;
	//cout << num << endl << m << endl << n << endl << ii << endl << jj << endl;

	return num;

}//////////////////////////////////////////////////

//*/

void slidewindow(int loop)
{
	 int width, highth;

	float imgcut[200][200] = { 0 };
	float feature[10000];// Feature[1000][10000];/////////////////////
	float *pp, *feat;
	float img[500][500] = { 0 };
	int winnum;
	 int totalnum;

	int slidepix = 5;
	//int width = 70;
	//int highth = 25;
	 int M;
	 int N;
         float part1schp[4];
 float part2schp[4];
 float part3schp[4];
 int partloop;
	int	windowl = width * 2 + 1;
	int windowh = highth * 2 + 1;
	int num, imgnum;
	int ii, jj;

	int step = 16;
	int oriy = 0, orix = 0;
	string picname;
	stringstream ss;
	ss << loop;
	ss >> picname;


	cout << "picname" << picname << endl;

	int m = windowh;
	int n = windowl;
	float moveinf;
	float mm, nn;
//char pic='1';
	ifstream inf;

    char imgdata[20];
    char str[20]=".txt";
char picdocmt[25]="pic/";
sprintf(imgdata,"%d",loop) ;
strcat(picdocmt,imgdata);  
strcat(picdocmt,str); 
  cout<<picdocmt<<endl;
inf.open(picdocmt);

	inf >> mm;
	inf >> nn;

        inf>>part1schp[0];
inf>>part1schp[1];
inf>>part1schp[2];
inf>>part1schp[3];

inf>>part2schp[0];
inf>>part2schp[1];
inf>>part2schp[2];
inf>>part2schp[3];

inf>>part3schp[0];
inf>>part3schp[1];
inf>>part3schp[2];
inf>>part3schp[3];
	M = int(mm);
	N = int(nn);
int globlesch=0;
	cout << "M" << M << "N" << N << endl;
	for (int i = 0; i<N - 14; i++)
	{
		inf >> moveinf;
		//cout<<moveinf<<endl;
	}

if (globlesch==1)
{
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			inf >> img[i][j];


		}
		//cout << endl;
	}
}


if (globlesch==0)
{

if (partloop==1)
{
for (int i = 0; i < part1schp[1]-highth-1; i++)
{
for (int j = 0; j<N; j++)
{
inf >> moveinf;
}
}



for (int i = 0; i < part1schp[0]-width-1; i++)
	{inf >> moveinf;
}
for (int i = 0; i < part1schp[3]-part1schp[1]+1+2*highth; i++)
{
for (int j = 0; j < part1schp[2]-part1schp[0]+1+2*width; j++)
	{inf >> img[i][j];
}

for (int k = 0; k < N-(part1schp[2]-part1schp[0]+1+2*width); k++)
	{inf >> moveinf;
}
}
M=part1schp[3]-part1schp[1]+1+2*highth;
N=part1schp[2]-part1schp[0]+1+2*width;

}

if (partloop==2)
{
for (int i = 0; i < part2schp[1]-highth-1; i++)
{
for (int j = 0; j<N; j++)
{
inf >> moveinf;
}
}


for (int i = 0; i < part2schp[0]-width-1; i++)
	{inf >> moveinf;
}
for (int i = 0; i < part2schp[3]-part2schp[1]+1+2*highth; i++)
{
for (int j = 0; j < part2schp[2]-part2schp[0]+1+2*width; j++)
	{inf >> img[i][j];
}

for (int k = 0; k < N-(part2schp[2]-part2schp[0]+1+2*width); k++)
	{inf >> moveinf;
}
}
M=part2schp[3]-part2schp[1]+1+2*highth;
N=part2schp[2]-part2schp[0]+1+2*width;

}

if (partloop==3)
{
for (int i = 0; i < part3schp[1]-highth-1; i++)
{
for (int j = 0; j<N; j++)
{
inf >> moveinf;
}
}


for (int i = 0; i < part3schp[0]-width-1; i++)
	{inf >> moveinf;
}
for (int i = 0; i < part3schp[3]-part3schp[1]+1+2*highth; i++)
{
for (int j = 0; j < part3schp[2]-part3schp[0]+1+2*width; j++)
	{inf >> img[i][j];
}

for (int k = 0; k < N-(part3schp[2]-part3schp[0]+1+2*width); k++)
	{inf >> moveinf;
}
}
M=part3schp[3]-part3schp[1]+1+2*highth;
N=part3schp[2]-part3schp[0]+1+2*width;

}


}
cout<<"M:"<<M<<endl<<"N:"<<N<<endl;

	inf.close();

	//cout << "img" << img[0][0] << endl << img[M - 1][N - 1] << endl;




	//计算样本量
	for (ii = 0; ii < ((M - windowh) / slidepix); ii++)
	{

		orix = 0;
		for (jj = 0; jj < ((N - windowl) / slidepix); jj++)
		{
			orix = orix + slidepix;


		}
		oriy = oriy + slidepix;
	}

	totalnum = ii*jj;
	cout << "totalnum" << totalnum << endl;



	//计算窗口维度
	ii = 0;
	jj = 0;
	for (int i = 1; i <= windowh - step; i = i + step)
	{
		//如果处理的m/step不是整数，最好是i=1:step:m-step
		jj = 0;
		for (int j = 1; j <= windowl - step; j = j + step)      //注释同上
		{
			jj = jj + 1;

		}
		ii = ii + 1;

	}
	//winnum = (ii - 1)*(jj - 1) * 36;
        winnum = (ii)*(jj) * 9;
cout<<"ii:"<<ii<<"jj:"<<jj<<endl;
	cout << "winnum" << winnum << endl;


	ofstream fout;
	fout.open("SRF/53855640/Featureout.data");
	fout << "double" << endl << totalnum << " " << winnum << endl << "dense" << endl;
	feat = &feature[0];
	//cout << img[0][8] << endl;
	//hogfun1(m, n,pp);

	//cout << (M - windowh) / slidepix << endl;
	oriy = 0;
	orix = 0;

	for (ii = 0; ii < ((M - windowh) / slidepix); ii++)
	{

		orix = 0;
		for (jj = 0; jj < ((N - windowl) / slidepix); jj++)
		{
			for (int p = 0; p < windowh; p++)
			{


				for (int q = 0; q < windowl; q++)
				{
					imgcut[p][q] = img[oriy + p][orix + q];
				}
			}
			pp = &imgcut[0][0];
			//cout << "test" << endl<<*pp<< endl<<*feat<<endl<<m<<endl<<n<<endl;
			num = hogfun1(m, n, pp, feat);

			for (int t = 0; t < num; t++)
			{

				//Feature[(ii - 1)*((N - windowl) / slidepix) + jj][t] = feature[t];
				fout << feature[t] << " ";
				//cout << feature[t] << endl;
			}
			fout << endl;


			orix = orix + slidepix;


		}
		oriy = oriy + slidepix;
	}

	imgnum = ii*jj;

	/*for (int in = 0; in < imgnum; in++)
	{
	for (int jn = 0; jn < num; jn++)
	{
	fout << Feature[in][jn];
	}
	fout << endl;

	}*/
	fout.close();

	//cout << "imgnum" << imgnum << endl<<"num"<<num<<endl;
	int label = 1;
	ofstream foutlabel;
	foutlabel.open("SRF/53855640/Featureoutlabel.label");
	foutlabel << "int" << endl << totalnum << " " << label << endl << "dense" << endl;
	for (int i = 0; i < totalnum; i++)
	{
		foutlabel << label - 1 << endl;
	}
	foutlabel.close();

}//////////////////////////////////////////////////////////////////
void printUsage()
{
	cout << "Usage: randomForest <config.conf> <option>" << endl;
	cout << "Options: " << "-train -test -trainAndTest" << endl;
}

void printTestUsage()
{
	cout << "Usage: randomForest <config.conf> -test <classifier.xml>" << endl;
}

HyperParameters parseHyperParameters(const std::string& confFile)
{
	HyperParameters hp;
	Config configFile;

	configFile.readFile(confFile.c_str());

	// DATA
	hp.trainData = (const char*)configFile.lookup("Data.trainData");
	hp.trainLabels = (const char*)configFile.lookup("Data.trainLabels");

hp.trainData1 = (const char*)configFile.lookup("Data.trainData1");
	hp.trainLabels1 = (const char*)configFile.lookup("Data.trainLabels1");

hp.trainData2 = (const char*)configFile.lookup("Data.trainData2");
	hp.trainLabels2 = (const char*)configFile.lookup("Data.trainLabels2");

hp.trainData3 = (const char*)configFile.lookup("Data.trainData3");
	hp.trainLabels3 = (const char*)configFile.lookup("Data.trainLabels3");
	hp.testData = (const char*)configFile.lookup("Data.testData");

	hp.testLabels = (const char*)configFile.lookup("Data.testLabels");
        

	hp.numLabeled = configFile.lookup("Data.numLabeled");
	hp.numClasses = configFile.lookup("Data.numClasses");
	hp.numTestNum = configFile.lookup("Data.numTestNum");
	hp.numTestNum1 = configFile.lookup("Data.numTestNum1");
	hp.numTestNum2 = configFile.lookup("Data.numTestNum2");
	hp.numTestNum3 = configFile.lookup("Data.numTestNum3");
        //hp.numTestNum4 = configFile.lookup("Data.numTestNum4");
	hp.weidu = configFile.lookup("Data.weidu");
	hp.width = configFile.lookup("Data.width");
	hp.highth = configFile.lookup("Data.highth");
	hp.width1 = configFile.lookup("Data.width1");
	hp.highth1 = configFile.lookup("Data.highth1");
	hp.width2 = configFile.lookup("Data.width2");
	hp.highth2 = configFile.lookup("Data.highth2");
	hp.width3 = configFile.lookup("Data.width3");
	hp.highth3 = configFile.lookup("Data.highth3");
         hp.width4 = configFile.lookup("Data.width4");
	hp.highth4 = configFile.lookup("Data.highth4");
hp.width5 = configFile.lookup("Data.width5");
	hp.highth5 = configFile.lookup("Data.highth5");
hp.width6 = configFile.lookup("Data.width6");
	hp.highth6 = configFile.lookup("Data.highth6");
hp.width7 = configFile.lookup("Data.width7");
	hp.highth7 = configFile.lookup("Data.highth7");
hp.width8 = configFile.lookup("Data.width8");
	hp.highth8 = configFile.lookup("Data.highth8");
hp.width9 = configFile.lookup("Data.width9");
	hp.highth9 = configFile.lookup("Data.highth9");
hp.width10 = configFile.lookup("Data.width10");
	hp.highth10 = configFile.lookup("Data.highth10");

	hp.winm = configFile.lookup("Data.winm");
	hp.winn = configFile.lookup("Data.winn");
	hp.picnum = configFile.lookup("Data.picnum");
        hp.picnum0= configFile.lookup("Data.picnum0");
hp.partnum= configFile.lookup("Data.partnum");
	/////////////////////////////////////////
	hp.className = (const char*)configFile.lookup("Data.className");
	hp.className1 = (const char*)configFile.lookup("Data.className1");
	hp.className2 = (const char*)configFile.lookup("Data.className2");
	hp.className3 = (const char*)configFile.lookup("Data.className3");
        hp.className4 = (const char*)configFile.lookup("Data.className4");
hp.className5 = (const char*)configFile.lookup("Data.className5");
hp.className6 = (const char*)configFile.lookup("Data.className6");
hp.className7 = (const char*)configFile.lookup("Data.className7");
hp.className8 = (const char*)configFile.lookup("Data.className8");
hp.className9 = (const char*)configFile.lookup("Data.className9");
hp.className10 = (const char*)configFile.lookup("Data.className10");
	// TREE
	hp.maxTreeDepth = configFile.lookup("Tree.maxDepth");
	hp.bagRatio = configFile.lookup("Tree.bagRatio");
	hp.numRandomFeatures = configFile.lookup("Tree.numRandomFeatures");
	hp.numProjFeatures = configFile.lookup("Tree.numProjFeatures");
	hp.useRandProj = configFile.lookup("Tree.useRandProj");
	hp.useGPU = configFile.lookup("Tree.useGPU");
	hp.useSubSamplingWithReplacement = configFile.lookup("Tree.subSampleWR");
	hp.verbose = configFile.lookup("Tree.verbose");
	hp.useInfoGain = configFile.lookup("Tree.useInfoGain");


	// FOREST
	hp.numTrees = configFile.lookup("Forest.numTrees");
	hp.useSoftVoting = configFile.lookup("Forest.useSoftVoting");
	hp.saveForest = configFile.lookup("Forest.saveForest");

	// OUTPUT
	hp.saveName = (const char *)configFile.lookup("Output.saveName");
	hp.savePath = (const char *)configFile.lookup("Output.savePath");
	hp.loadName = (const char *)configFile.lookup("Output.loadName");

	return hp;
}

void saveForest(const Forest& forest, const std::string& filename)
{
	/*
	const xmlNodePtr rootNode = forest.save();
	xmlDocPtr doc = xmlNewDoc( reinterpret_cast<const xmlChar*>( "1.0" ) );
	xmlDocSetRootElement( doc, rootNode );
	xmlSaveFormatFileEnc( filename.c_str(), doc, "UTF-8", 1 );
	xmlFreeDoc( doc );
	*/
	FILE *fp = fopen(filename.c_str(), "w");
	fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
	forest.save(fp);
}
float Ix[200][200], Ied[200][200], Iphase[200][200];
float Iy[200][200];// Ied[500][500], Iphase[500][500];
float img[200][200];
float Cell[9][1000][1000];
float tmpx[26][26], tmped[26][26], tmpphase[26][26], f[9], sum1 = 0;
float imginput[110][110];
int width, highth;
int partloop;
float part1schp[4],part2schp[4],part3schp[4];
int totalnum;
float evaconf[200];
float evaconf1[200];
int M, N;
/*
int main()
{  
 return 0;
}
*/

int main(int argc, char** argv)
{  
 float evaconf[200];

  for(int i=0;i<200;i++)
{
      evaconf[i]=0;
}
	// at least, we need two arguments
	if (argc < 2) {
		cout << "----------------------------------------------------------------------------------------------------" << endl;
		cout << "Illegal arguments specified. Your options:" << endl;
		cout << "----------------------------------------------------------------------------------------------------" << endl;
		printUsage();
		printTestUsage();
		cout << "----------------------------------------------------------------------------------------------------" << endl;
		return -1;
	}

	// Train Random Forest
	if (argc == 3 && (std::string(argv[2]) == "-train"))
	{
		const std::string confFile(argv[1]);
		HyperParameters hp = parseHyperParameters(confFile);
		Forest rf(hp);
		//Tree t(hp);

		FileData fileData(hp.trainData, hp.trainLabels);
		fileData.readData();
		fileData.readLabels();

		rf.train(fileData.getData(), fileData.getLabels());
		//saveForest(rf,"classifier.xml");
		saveForest(rf, hp.className);
	}
         else if (std::string(argv[2]) == "-trainallpart")
        {
                const std::string confFile(argv[1]);
		HyperParameters hp = parseHyperParameters(confFile);
		Forest rf1(hp);
		//Tree t(hp);
Forest rf2(hp);
Forest rf3(hp);
		FileData fileData1(hp.trainData1, hp.trainLabels1);
		fileData1.readData();
		fileData1.readLabels();

		rf1.train(fileData1.getData(), fileData1.getLabels());
		//saveForest(rf,"classifier.xml");
		saveForest(rf1, hp.className1);

                FileData fileData2(hp.trainData2, hp.trainLabels2);
		fileData2.readData();
		fileData2.readLabels();

		rf2.train(fileData2.getData(), fileData2.getLabels());
		//saveForest(rf,"classifier.xml");
		saveForest(rf2, hp.className2);


                FileData fileData3(hp.trainData3, hp.trainLabels3);
		fileData3.readData();
		fileData3.readLabels();

		rf3.train(fileData3.getData(), fileData3.getLabels());
		//saveForest(rf,"classifier.xml");
		saveForest(rf3, hp.className3);

               // FileData fileData3(hp.trainData3, hp.trainLabels3);
		//fileData3.readData();
		//fileData3.readLabels();

		//rf3.train(fileData3.getData(), fileData3.getLabels());
		//saveForest(rf,"classifier.xml");
		//saveForest(rf3, hp.className3);

        }
	// Test Random Forest
	else if (std::string(argv[2]) == "-test")
	{


		if (argc == 3)
		{
			const std::string confFile(argv[1]);
			HyperParameters hp = parseHyperParameters(confFile);
			cout << "hp" << endl;////////////////////

			width = hp.width;
			highth = hp.highth;
			//M=hp.winm;
			//N=hp.winn;
			int loop = 3;
			slidewindow(loop);
			Forest rf(hp, hp.className);
			cout << "rf" << endl;////////////////////
                        
			//cout<<"argv"<<argv[3]<<endl;
			//cout<<"hp"<<hp.className<<endl;
			FileData testData(hp.testData, hp.testLabels);
			testData.readData();
			testData.readLabels();
			cout << "readdata" << endl;////////////////////
			rf.eval(testData.getData(), testData.getLabels());

			///////////////////////

			rf.evaConf(rf.m_confidences, totalnum, hp.numClasses, width, highth, M, N,1);
		}
		else
		{
			printTestUsage();
		}
	}
	// Test Random Forest with different part
	

	else if (std::string(argv[2]) == "-testallpart")
	{
		if (argc == 3)
		{       float partmaxeva;
                        int partmaxnum;
 int partloop;
int part1errornum=0;
int part2errornum=0;
int part3errornum=0;
int part4errornum=0;
int parterrornum[30];

			//float evaconf1[200];
			//float evaconf2[200];
			//float evaconf3[200];
                        //float evaconf4[200];
float evaconfpart[30][200]={0};
float evaconf1part[30][200]={0};
char classnm[30]="hp.className";
char str[10];
int prtnum=1;
sprintf(str,"%d",prtnum);
strcat(classnm,str);
cout<<"classnm:"<<classnm<<endl;
			const std::string confFile(argv[1]);
cout<<"test"<<endl;
			HyperParameters hp = parseHyperParameters(confFile);

			Forest rf1(hp, hp.className1);

			Forest rf2(hp, hp.className2);

			Forest rf3(hp, hp.className3);

                        Forest rf4(hp, hp.className4);

Forest rf5(hp, hp.className5);

			Forest rf6(hp, hp.className6);

			Forest rf7(hp, hp.className7);

                        Forest rf8(hp, hp.className8);

Forest rf9(hp, hp.className9);

			Forest rf10(hp, hp.className10);

			
ofstream foutfinal;
	foutfinal.open("SRF/53855640/foutfinal.txt");
int finalnum[hp.picnum-hp.picnum0+1];
int carplus1=1;
int carnum1=1;


			for (int loop = hp.picnum0; loop <= hp.picnum; loop++)
			{foutfinal<<endl<<endl<<loop<<endl<<endl<<endl;
////////////////////////////part1///////////////////////////
for (partloop=1;partloop<=hp.partnum;partloop++)
{
if (partloop==1)
{				width = hp.width1;
				highth = hp.highth1;
}
if (partloop==2)
{				width = hp.width2;
				highth = hp.highth2;
}
if (partloop==3)
{				width = hp.width3;
				highth = hp.highth3;
}
if (partloop==4)
{				width = hp.width4;
				highth = hp.highth4;
}
if (partloop==5)
{				width = hp.width5;
				highth = hp.highth5;
}
if (partloop==6)
{				width = hp.width6;
				highth = hp.highth6;
}
if (partloop==7)
{				width = hp.width7;
				highth = hp.highth7;
}
if (partloop==8)
{				width = hp.width8;
				highth = hp.highth8;
}
if (partloop==9)
{				width = hp.width9;
				highth = hp.highth9;
}
if (partloop==10)
{				width = hp.width10;
				highth = hp.highth10;
}
				//M=hp.winm;
				//N=hp.winn;
				slidewindow(loop);



				FileData testData(hp.testData, hp.testLabels);
				testData.readData();
				testData.readLabels();
if (partloop==1)
{			        rf1.getConfidences();
				rf1.eval(testData.getData(), testData.getLabels());


				rf1.evaConf(rf1.m_confidences, totalnum, hp.numClasses, width, highth, M, N,1);

}
if (partloop==2)
{			        rf2.getConfidences();
				rf2.eval(testData.getData(), testData.getLabels());


				rf2.evaConf(rf2.m_confidences, totalnum, hp.numClasses, width, highth, M, N,2);

}
if (partloop==3)
{			        rf3.getConfidences();
				rf3.eval(testData.getData(), testData.getLabels());


				rf3.evaConf(rf3.m_confidences, totalnum, hp.numClasses, width, highth, M, N,3);

}
if (partloop==4)
{			        rf4.getConfidences();
				rf4.eval(testData.getData(), testData.getLabels());


				rf4.evaConf(rf4.m_confidences, totalnum, hp.numClasses, width, highth, M, N,4);

}

if (partloop==5)
{			        rf5.getConfidences();
				rf5.eval(testData.getData(), testData.getLabels());


				rf5.evaConf(rf5.m_confidences, totalnum, hp.numClasses, width, highth, M, N,5);

}
if (partloop==6)
{			        rf6.getConfidences();
				rf6.eval(testData.getData(), testData.getLabels());


				rf6.evaConf(rf6.m_confidences, totalnum, hp.numClasses, width, highth, M, N,6);

}
if (partloop==7)
{			        rf7.getConfidences();
				rf7.eval(testData.getData(), testData.getLabels());


				rf7.evaConf(rf7.m_confidences, totalnum, hp.numClasses, width, highth, M, N,7);

}
if (partloop==8)
{			        rf8.getConfidences();
				rf8.eval(testData.getData(), testData.getLabels());


				rf8.evaConf(rf8.m_confidences, totalnum, hp.numClasses, width, highth, M, N,8);

}
if (partloop==9)
{			        rf9.getConfidences();
				rf9.eval(testData.getData(), testData.getLabels());


				rf9.evaConf(rf9.m_confidences, totalnum, hp.numClasses, width, highth, M, N,9);

}
if (partloop==10)
{			        rf10.getConfidences();
				rf10.eval(testData.getData(), testData.getLabels());


				rf10.evaConf(rf10.m_confidences, totalnum, hp.numClasses, width, highth, M, N,10);

}
				for (int i = 0; i < 200; i++)
				{	evaconfpart[partloop-1][i] = evaconf[i];

evaconf1part[partloop-1][i] = evaconf1[i];
if (i<6)
cout<<evaconf1part[partloop-1][i]<<endl;

}


}








				//  cout<<"eva"<<evaconf3[4]<<endl;
////////////////////////////compute///////////////////////////
				float evaconf[200];

				int maxnum = 1;
				float max;
				cout << "totalresult" << endl;

				for (int i = 0; i<hp.numClasses-1; i++)
				{
//					evaconf[i] = (evaconfpart[0][i] + evaconfpart[1][i] + evaconfpart[2][i]+evaconfpart[3][i]+evaconfpart[4][i]+evaconfpart[5][i]+evaconfpart[6][i]+evaconfpart[7][i]+evaconfpart[8][i]+evaconfpart[9][i]) / hp.partnum;


evaconf[i] = (evaconf1part[0][i] + evaconf1part[1][i] + evaconf1part[2][i]+evaconf1part[3][i]+evaconf1part[4][i]+evaconf1part[5][i]+evaconf1part[6][i]+evaconf1part[7][i]+evaconf1part[8][i]+evaconf1part[9][i]) / hp.partnum;
					if (i == 0)
					{
						max = evaconf[0];
                                           if(evaconf[0]!=0)
						maxnum = 1;
                                           if(evaconf[0]==0)
						maxnum = 0;
					}
					else if (evaconf[i]>max)
					{
						max = evaconf[i];
						maxnum = i + 1;
					}
					cout << i + 1 << " " << evaconf[i] << endl;

				}

finalnum[loop-hp.picnum0]=maxnum;

				cout << "Car Classification Evaculation Result:" << maxnum << " " << "confidence:" << max << endl;
				cout << "totalresult" << endl;
				for (int i = 0; i < hp.numClasses-1; i++)
				{
					//evaconf[i] = (evaconf1[i] + evaconf2[i] + evaconf3[i] +evaconf4[i])/ 4;
					cout << evaconf[i] << endl;
foutfinal<<i+1<<"  "<<evaconf[i]<<endl;
				}

foutfinal<<endl<<endl<<"evaclass  "<<maxnum<<endl<<endl<<endl;


////////////////////////////parterrornum///////////////////////////

char str[20];
char part[20]="part";
char parterrornumm[20]="parterrornum";
char oriparterrornumm[20]="parterrornum";


for(int partloop=1;partloop<=hp.partnum;partloop++)
{strcpy(parterrornumm,oriparterrornumm);

sprintf(str,"%d",partloop);
strcat(part,str);
strcat(parterrornumm,str);
foutfinal<<part<<endl<<endl;

                          for (int i = 0; i < hp.numClasses-1; i++)
				{if(i==0)
                                 {partmaxeva=evaconfpart[partloop-1][i];
                                  partmaxnum=i;}
				if (evaconfpart[partloop-1][i]>partmaxeva)
                                   {partmaxeva=evaconfpart[partloop-1][i];
                                  partmaxnum=i;}
					
foutfinal<<i+1<<"  "<<evaconfpart[partloop-1][i]<<endl;
				}
foutfinal<<parterrornumm<<" "<<partmaxnum+1<<endl;




if (carplus1==6)
{
carnum1=carnum1+1;
carplus1=1;
}
if((partmaxnum+1)!=carnum1)
{parterrornum[partloop-1]=parterrornum[partloop-1]+1;

}



}


carplus1=carplus1+1;

			}
foutfinal<<" finalresult: "<<endl;

int errornum=0;
int carplus=1;
int carnum=1;
for (int i = 0; i < hp.picnum-hp.picnum0+1; i++)
				{
					
					


//if(finalnum[i]!=hp.picnum0+i)
//errornum++;







if (carplus==6)
{
carnum=carnum+1;
carplus=1;
}

if(finalnum[i]!=carnum)
errornum=errornum+1;


carplus=carplus+1;


//foutfinal<<i+hp.picnum0<<"  "<<finalnum[i]<<endl;

foutfinal<<i+hp.picnum0<<" "<<carnum<<"  "<<finalnum[i]<<endl;



				}

float errorrate=float(errornum)/float(hp.picnum-hp.picnum0+1);
foutfinal<<"errornum:"<<errornum<<endl<<"errorrate:"<<errorrate<<endl<<"correctrate:"<<(1-errorrate)*100<<"%"<<endl;



/////////////////////////computeparterrornum errorrate//////////////////////////


char parterrornumm1[20]="parterrornum";
char oriparterrornumm1[20]="parterrornum";

char parterrorrate[20]="partcorrectrate";
char oriparterrorrate[20]="partcorrectrate";

for(int partloop=1;partloop<=hp.partnum;partloop++)
{strcpy(parterrornumm1,oriparterrornumm1);
strcpy(parterrorrate,oriparterrorrate);
sprintf(str,"%d",partloop);

strcat(parterrornumm1,str);
strcat(parterrorrate,str);



foutfinal<<parterrorrate<<":"<<(1-float(parterrornum[partloop-1])/float(hp.picnum-hp.picnum0+1))*100<<"%"<<endl;


foutfinal<<parterrornumm1<<":"<<parterrornum[partloop-1]<<endl;
}

foutfinal.close();

		}

		else
		{
			printTestUsage();
		}
	}


else if (std::string(argv[2]) == "-testback")
	{

//string hello("hello,world!");
//int arr[]={1,2,4,4};
//BOOST_FOREACH(int i,arr)
//{

//    cout<<i<<endl;


//}









		if (argc == 3)
		{       float partmaxeva;
                        int partmaxnum;
int part1errornum=0;

matrix<float> test_confidences;
			float evaconf1[200];
			
			const std::string confFile(argv[1]);
			HyperParameters hp = parseHyperParameters(confFile);
			Forest rf1(hp, hp.className);//xmlFreeDoc( forestDoc );
		cout<<"loadover"<<endl;
ofstream foutfinal;
	foutfinal.open("SRF/53855640/foutfinal.txt");
int finalnum[hp.picnum];
			
				width = hp.width1;
				highth = hp.highth1;
				//M=hp.winm;
				//N=hp.winn;
				//slidewindow(loop);



				FileData testData(hp.trainData, hp.trainLabels);
				testData.readData();
				testData.readLabels();
				rf1.getConfidences();//return my_confidences
//////////////////////////////////step1/////////////////////////////////
cout<<"step1"<<endl;

				rf1.eval(testData.getData(), testData.getLabels());

test_confidences=rf1.getConfidences();
//scout<<"testconfidences"<<test_confidences(1,1)<<endl;
				rf1.evaConf(rf1.m_confidences, hp.numLabeled, hp.numClasses, width, highth, M, N,1);
				




ofstream fout;
    fout.open("myconfidence.txt");


ifstream fin;
	fin.open("testback.label");



int maxx[200000];
float maxconf;
int temp;
float errornum=0;
float errorrate;
for (int i = 0; i < hp.numLabeled; i++)
{foutfinal<<i+1<<" ";
    maxx[i]=0;
maxconf=0;
for (int j = 0; j < hp.numClasses; j++)
{if(rf1.m_confidences(i,j)>maxconf)
{maxx[i]=j;
maxconf=rf1.m_confidences(i,j);
}
   fout<<rf1.m_confidences(i,j)<<" ";
}

fout<<endl;
fin>>temp;
foutfinal<<maxx[i]<<" "<<temp<<" "<<maxconf<<" "<<rf1.m_confidences(i,temp)<<endl;
if(maxx[i]!=temp)
{
errornum=errornum+1;

}

}



errorrate=errornum/float(hp.numLabeled);
foutfinal<<"errorrate"<<errorrate<<endl<<"errornum"<<errornum<<endl;
fin.close();
foutfinal.close();
fout.close();


		}

		else
		{
			printTestUsage();
		}
	}
else if (argc == 3 && (std::string(argv[2]) == "-trainandtestallpart"))

{
    ////////////////////////////train
 const std::string confFile(argv[1]);
		HyperParameters hp = parseHyperParameters(confFile);
		Forest rf1(hp);
		//Tree t(hp);
Forest rf2(hp);
Forest rf3(hp);
		FileData fileData1(hp.trainData1, hp.trainLabels1);
		fileData1.readData();
		fileData1.readLabels();

		rf1.train(fileData1.getData(), fileData1.getLabels());
		//saveForest(rf,"classifier.xml");
		saveForest(rf1, hp.className1);

                FileData fileData2(hp.trainData2, hp.trainLabels2);
		fileData2.readData();
		fileData2.readLabels();

		rf2.train(fileData2.getData(), fileData2.getLabels());
		//saveForest(rf,"classifier.xml");
		saveForest(rf2, hp.className2);

                FileData fileData3(hp.trainData3, hp.trainLabels3);
		fileData3.readData();
		fileData3.readLabels();

		rf3.train(fileData3.getData(), fileData3.getLabels());
		//saveForest(rf,"classifier.xml");
		saveForest(rf3, hp.className3);

////////////////////////train

///////////////////test
if (argc == 3)
		{
			float evaconf1[200];
			float evaconf2[200];
			float evaconf3[200];
			//const std::string confFile(argv[1]);
			//HyperParameters hp = parseHyperParameters(confFile);
			Forest rf1(hp, hp.className1);
			Forest rf2(hp, hp.className2);
			Forest rf3(hp, hp.className3);
ofstream foutfinal;
	foutfinal.open("SRF/53855640/foutfinal.txt");
int finalnum[hp.picnum];
			for (int loop = 1; loop <= hp.picnum; loop++)
			{foutfinal<<endl<<endl<<loop<<endl<<endl<<endl;
				width = hp.width1;
				highth = hp.highth1;
				//M=hp.winm;
				//N=hp.winn;
				slidewindow(loop);



				FileData testData(hp.testData, hp.testLabels);
				testData.readData();
				testData.readLabels();
				rf1.getConfidences();
				rf1.eval(testData.getData(), testData.getLabels());


				rf1.evaConf(rf1.m_confidences, totalnum, hp.numClasses, width, highth, M, N,1);
				for (int i = 0; i < 200; i++)
					evaconf1[i] = evaconf[i];

				width = hp.width2;
				highth = hp.highth2;

				slidewindow(loop);

				testData.readData();
				testData.readLabels();
				rf2.getConfidences();
				rf2.eval(testData.getData(), testData.getLabels());
				rf2.evaConf(rf2.m_confidences, totalnum, hp.numClasses, width, highth, M, N,2);
				for (int i = 0; i < 200; i++)
					evaconf2[i] = evaconf[i];

				width = hp.width3;
				highth = hp.highth3;
				cout << "wid" << width << "highth" << highth << endl;
				slidewindow(loop);

				testData.readData();
				testData.readLabels();
				rf3.getConfidences();
				rf3.eval(testData.getData(), testData.getLabels());
				//cout << "totalnum" << totalnum << endl;
				//cout << "M" << M << "N" << N << endl;
				rf3.evaConf(rf3.m_confidences, totalnum, hp.numClasses, width, highth, M, N,3);

				for (int i = 0; i < 200; i++)
					evaconf3[i] = evaconf[i];
				//  cout<<"eva"<<evaconf3[4]<<endl;
				////////////////compute
				float evaconf[200];
				int maxnum = 1;
				float max;
				cout << "totalresult" << endl;

				for (int i = 0; i<hp.numClasses-1; i++)
				{
					evaconf[i] = (evaconf1[i] + evaconf2[i] + evaconf3[i]) / 3;
					if (i == 0)
					{
						max = evaconf[0];
                                           if(evaconf[0]!=0)
						maxnum = 1;
                                           if(evaconf[0]==0)
						maxnum = 0;
					}
					else if (evaconf[i]>max)
					{
						max = evaconf[i];
						maxnum = i + 1;
					}
					cout << i + 1 << " " << evaconf[i] << endl;

				}

finalnum[loop]=maxnum;

				cout << "Car Classification Evaculation Result:" << maxnum << " " << "confidence:" << max << endl;
				cout << "totalresult" << endl;
				for (int i = 0; i < hp.numClasses-1; i++)
				{
					evaconf[i] = (evaconf1[i] + evaconf2[i] + evaconf3[i]) / 3;
					cout << evaconf[i] << endl;
foutfinal<<i+1<<"  "<<evaconf[i]<<endl;
				}

foutfinal<<endl<<endl<<"evaclass  "<<maxnum<<endl<<endl<<endl;
foutfinal<<"part1"<<endl<<endl;
                          for (int i = 0; i < hp.numClasses-1; i++)
				{
					
					
foutfinal<<i+1<<"  "<<evaconf1[i]<<endl;
				}
foutfinal<<endl<<"part2"<<endl<<endl;
for (int i = 0; i < hp.numClasses-1; i++)
				{
					
					
foutfinal<<i+1<<"  "<<evaconf2[i]<<endl;
				}
foutfinal<<endl<<"part3"<<endl<<endl;
for (int i = 0; i < hp.numClasses-1; i++)
				{
					
					
foutfinal<<i+1<<"  "<<evaconf3[i]<<endl;
				}


			}
foutfinal<<" finalresult: "<<endl;

int errornum=0;
for (int i = 0; i < hp.picnum; i++)
				{
					
					
foutfinal<<i+1<<"  "<<finalnum[i+1]<<endl;

if(finalnum[i+1]!=i+1)
errornum++;
				}

float errorrate=float(errornum)/float(hp.picnum);
foutfinal<<"errornum:"<<errornum<<endl<<"errorrate:"<<errorrate<<endl;
foutfinal.close();

		}

		else
		{
			printTestUsage();
		}
///////////////////////test



}













	// Train a Tree
	else if (argc == 3 && (std::string(argv[2]) == "-trainAndTestTree"))
	{
		const std::string confFile(argv[1]);
		HyperParameters hp = parseHyperParameters(confFile);

		Tree rt(hp);

		timeIt(1);
		FileData trainData(hp.trainData, hp.trainLabels);
		trainData.readData();
		trainData.readLabels();

		rt.train(trainData.getData(), trainData.getLabels());
		cout << "\tTraining time = " << timeIt(0) << " seconds" << endl;

		std::vector<int> nodeLabels = rt.getNodeLabels();
		cout << "Classes which has at least a leaf node ... " << endl;
		dispVector(nodeLabels);

		timeIt(1);
		FileData testData(hp.testData, hp.testLabels);
		testData.readData();
		testData.readLabels();
		rt.eval(testData.getData(), testData.getLabels());
		cout << "\tTest time = " << timeIt(0) << " seconds" << endl;
	}
	// Train and Test Random Forest
	else if (argc == 3 && (std::string(argv[2]) == "-trainAndTest"))
	{
		const std::string confFile(argv[1]);
		HyperParameters hp = parseHyperParameters(confFile);

		Forest rf(hp);

		timeIt(1);
		FileData trainData(hp.trainData, hp.trainLabels);
		trainData.readData();
		trainData.readLabels();



		rf.train(trainData.getData(), trainData.getLabels());
		cout << "\tTraining time = " << timeIt(0) << " seconds" << endl;

		timeIt(1);
		FileData testData(hp.testData, hp.testLabels);
		testData.readData();
		testData.readLabels();
		rf.eval(testData.getData(), testData.getLabels());
		cout << "\tTest time = " << timeIt(0) << " seconds" << endl;
		if (hp.saveForest) {
			saveForest(rf, "testforest.xml");
		}
	}

	// Train and Test Random Forest
	else if (argc == 3 && (std::string(argv[2]) == "-trainPairsAndTest"))
	{
		const std::string confFile(argv[1]);
		HyperParameters hp = parseHyperParameters(confFile);

		PairForest rf(hp);

		timeIt(1);
		FileData trainData(hp.trainData, hp.trainLabels);
		trainData.readData();
		trainData.readLabels();
		trainData.createPairs();

		rf.train(trainData.getPairs());
		cout << "\tTraining time = " << timeIt(0) << " seconds" << endl;
		cout << "\tTrain result: ";
		rf.test(trainData.getPairs());


		timeIt(1);
		FileData testData(hp.testData, hp.testLabels);
		testData.readData();
		testData.readLabels();
		testData.createPairs();
		cout << "\tTest result: ";
		rf.test(testData.getPairs());

		cout << "\tTest time = " << timeIt(0) << " seconds" << endl;

	}


	return 0;
}

