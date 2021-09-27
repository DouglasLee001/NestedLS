#ifndef _BASIC_H
#define _BASIC_H
#include <climits>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <map>
#include "Array.h"
#include "MyHeap.h"

using namespace std;

#define Pop(stack) stack[--stack##_fill_pointer]
#define Push(item, stack) stack[stack##_fill_pointer++] = item

//Debug control
#define NDEBUG
#define NDEBUG_OUTPUT 
#define TREE_SCORE
#define NORIGIN_SCORE


bool moduleRemoveRedundant = false;

void printDebugMsg(string msg)
{
#ifdef DEBUG_OUTPUT
    cout << msg << endl;
#endif
}
void printDebugRemove(int v, int step, int time)
{
#ifdef DEBUG_OUTPUT
    cout << "Remove: " << v << " at " << step << " since last time: " << time << endl;
#endif
}
void printDebugAdd(int v, int step, int time)
{
#ifdef DEBUG_OUTPUT
    cout << "Add: " << v << " at " << step << " since last time: " << time << endl;
#endif
}

typedef long long llong;
typedef unsigned int uint;

struct Edge
{
    int v1;
    int v2;
};
typedef struct Edge2
{
    int u, v;
    int next;
            
} EdgeLib;

int *first;

EdgeLib *edge;

chrono::steady_clock::time_point start;

enum MethodType
{
    CutTree = 0,
    Tarjan,
};
MethodType curMethod = MethodType::CutTree;

enum ChooseMode
{
    ModeA = 0,
    ModeB,
    ModeC,
    ModeD,
};

ChooseMode currentMode = ChooseMode::ModeA;

int *subWeight;
int *pre_deci_step;
int *temp_pre_deci_step;
int add_step = 0;
bool isUseBestSolutionInfo = false;
bool isRecordBestSolution = false;

int floor0 = 100;
int ceilingTimes = 10;
int insTimes = 2;
int tarjanRemoveSize = 50;
uint seed = 1;
int cutoff_time;
int tarjanCntUB = 8;
int improveRate = 200;
double reduceRate = 0.7;
long maxTrunkStepAction = 1000;
long noImproveTrunkTryStep = 100;
int bestWeightInTurn;

int floor1 = 10000;
double weightAdjustCoefficient = 0.7;
double weight_para_aphle = 0.7;
long weightThresholdCoefficient = 300;
long noImproveTryStep = 10000;
double restartUpdateFrequencyThreshold = 30.0;
int tarjanLoopCoefficient = 1;
int tarjanCntLB = 1;

int instance0; //50*2
int gap0;      //50*10
int instance1;
int gap1;

bool isForceRestart = false;
bool isShutDownSearchStep = false;
const double improveThresholdCoefficient = 0.0005;
int searchStepCnt = 0;
const int searchStepCntThreshold = 10;

llong step = 1;
llong stepAction = 1;
llong NOimprovementstep = 0;
llong steppre;
int try_step;
int betweemNoimpr = 1; //floor0*5*betweenNoimpr:[1,8]
llong Noimprovementstepfocus = 0;
int choosedremove_v, choosedadd_v;
int v_num;
int e_num;

//Edge *edge;
//int *edge_weight;
int maxNeighborSize = 3;
int *frequency;

double *weight;
long *subscore;
bool rightAfternewlow = true;
Array *removedNodeNeighbor;
Array *redundantNodes;
Array *leafArray;

long *score;
long *score_in_d;
llong *time_stamp;
llong *time_stamp_tree;
int *candidate;
int *index_in_candidate;
int candidate_size;

int *v_weight;
int **v_edges;
long **v_adj;
int *v_degree;
int *v_fixed;
int *fixedSet;
int fixedNum = 0;
int *v_threshold;
int *isInS; 
int *S;     
int Snum = 0;
int *Stack;
int *SF;
int *child;
long *onlydominate;

int *dominated;
int *dominatedByDApos;
long *onlyDominateByDApos;
Array *greyPointArray;
Array *undomPointArray;
Array *candidateArray;
Array *candidateAposArray;
Array *candidateForDApos;
Array *notDominatedByDApos;
Array *greyPointArrayDApos;
int c_size = 0;
int *v_in_c;
int *last_v_in_c;
int maxDegreeNode = 1;
llong now_weight;

int best_c_size;
int *best_v_in_c;
int *lastTurn_v_in_c;
int *best_dominated;
double best_comp_time;
llong best_step;

//Trunk
int curTrunkWeight = 0;
int bestTrunkWeight = 0;
long *trunkTabuRemove;
long *trunkTabuAdd;
llong trunkStep = 1;
llong trunkStepAction = 1;

long *taburemove;
long *tabuadd;
int tabutenue = 5;

int *conf_change;
int *pre;
int *CV;
int *WV;
//Tarjan
int *dnf;
int *low;
int *isCut;
//int *recordisCut;
//int recordcutIndex;
int *RemovedPoint;
int *AddedPoint;
int *cutPointSet;
int *initIndex;

int cutIndex;
int ave_weight;
int connectedNum = 0;
int minUndom = 0;
int undomafteradd;
llong averagedegree = 0;

double totalweight;
double weightthreshold;
double currentWeight = 0;
double bestWeight;

int tarjanLoopCount = 1;

int *father;
int *childnum;

MyHeap *constructSolutionHeap;
MyHeap *constructTreeHeap;

//build_free.h
double TimeElapsed();
int BuildInstance(string);
void FreeMemory();

//construct_initscore.h
// void init_increase_dominate(int, int);
// int find(int);
// void join(int, int);
// int calCV(int);
// void lowerScore();
// void initSubScore();
// void joinV(int);
// void addToS(int);
// void updateS(int);
// int chooseMax();
// void addNodeInit(int);
// void ConstructByInitScore();
// void restartIncreaseDominate(int, int);
// void restartAdd(int);

//construct_solution.h
void resetScore();
void constructIncreaseDominate(int, int);
void constructAdd(int, bool);
void prepareForLocalSearch();
void ConstructSolution();
void ConstructSolutionToBestWeight();
int NewSolutionChooseVFromMethodA();
int NewSolutionChooseVFromMethodB();
int NewSolutionChooseVFromMethodC();
int NewSolutionChooseVFromMethodD();
double lastRepeatDegree = 0.0;
double getNewSolutionRepeatDegree();
void treeIncreaseDominate(int, int, bool);
void treeDecreaseDominate(int, bool);

//update.h
int ind = 0;
int root = 1;
int maxScore = 0;
int maxPoint = 0;
int Toroot = -1;
void ResetCandidate();
void UpdateBestSolution();
void cutPoint(int, int);
void cutPointNoRecur(int);
void cutPoint1(long);
bool judgeCut(int);
void updateWeight();
void updateSubscore();
void ConstructRestartSolution();

//checker.h
bool test_score();
bool CheckSolution();
bool CheckSolutionIsConnected();
bool CheckGraphIsConnected();

//fastds-pure.h
inline void Undom(int);
inline void Dom(int);
void increase_dominate(long, long);
void decrease_dominate(int);
bool Add(int);
bool Remove(int);
void addToDApos(int, bool);
void removeFromDApos(int, bool);
void MarkCut();
void MarkCuttree();
void MarkCutTreeTrunk();
void addWeight(int);
void minusWeight(int, int);
void addUpdate(int, bool);
void removeUpdate(int, bool);
int ChooseAddVsubscorefast(int cnt, int choice);
int ChooseAddVsubscorefastAspration();
int ChooseAddVForConstructTrunk(const int);
int ChooseRemoveVTopof();
int ChooseRemoveVTopofBMS(const int count, int choice);
int ChooseRemoveVFromArray(Array *, int);
int ChooseRemoveVFromDApos(const int);
void localSearchTrunk();
void localSearchFramework1();
void localSearchFramework2();
void Framework1CutTree();
void Framework1Tarjan();
void Framework2CutTree();
void Framework2CutTreeFocus();
void Framework2TarjanFocus();
void Framework2TarjanFocus2();
void Framework2TarjanScatter();
void updateRedundantV(int);
void RemoveRedundant();
void Restart();
void enter_ls();
#endif
