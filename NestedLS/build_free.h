#ifndef _BUILD_FREE_H
#define _BUILD_FREE_H

#include "basic.h"

double TimeElapsed()
{
    chrono::steady_clock::time_point finish = chrono::steady_clock::now();
    chrono::duration<double> duration = finish - start;
    return duration.count();
}

int BuildInstance(string filename)
{
    instance0 = floor0 * insTimes; //50*2
    gap0 = floor0 * ceilingTimes;  //50*10
    instance1 = floor1 * insTimes;
    gap1 = floor1 * ceilingTimes;
    ifstream infile;
    infile.open(filename);
    if (infile.fail())
    {
        cout << "# Read file error: there is no such file existed!" << endl;
        getchar();
        exit(0);
    }
    char tempStr[100];
    infile >> tempStr;
    if (infile.eof())
    {
        cout << "# Read file error: file is empty!" << endl;
        getchar();
        exit(0);
    }

    int edgeCount = 0;
    while (!infile.eof())
    {
        if (strcmp(tempStr, "p") == 0)
        {
            infile >> tempStr;
            infile >> v_num >> e_num;
            /*****************new******************/
            frequency = new int[v_num + 1];
            weight = new double[v_num + 1];
            subWeight = new int[v_num + 1];
            pre_deci_step = new int[v_num + 1];
            candidate = new int[v_num + 1];
            index_in_candidate = new int[v_num + 1];
            temp_pre_deci_step = new int[v_num + 1];
            weightthreshold = (v_num + 1) * weightThresholdCoefficient;
            v_threshold = new int[v_num + 1];
            score = new long[v_num + 1];
            score_in_d = new long[v_num + 1];
            subscore = new long[v_num + 1];
            time_stamp = new llong[v_num + 1];
            time_stamp_tree = new llong[v_num + 1];
            greyPointArray = new Array(v_num);
            undomPointArray = new Array(v_num);
            candidateArray = new Array(v_num);
            leafArray = new Array(v_num);
            candidateAposArray = new Array(v_num);
            greyPointArrayDApos = new Array(v_num);
            first = new int[v_num + 1];
            v_edges = new int *[v_num + 1];
            v_adj = new long *[v_num + 1];
            v_degree = new int[v_num + 1];
            fixedSet = new int[v_num + 1];
            initIndex = new int[v_num + 1];
            v_in_c = new int[v_num + 1];
            last_v_in_c = new int[v_num + 1];
            lastTurn_v_in_c = new int[v_num + 1];
            best_v_in_c = new int[v_num + 1];
            conf_change = new int[v_num + 1];
            dominated = new int[v_num + 1];
            dominatedByDApos = new int[v_num + 1];
            best_dominated = new int[v_num + 1];
            isInS = new int[v_num + 1];
            edge = new EdgeLib[e_num + e_num];
            S = new int[v_num + 1];
            CV = new int[v_num + 1];
            WV = new int[v_num + 1];
            pre = new int[v_num + 1];
            dnf = new int[v_num + 1];
            low = new int[v_num + 1];
            isCut = new int[v_num + 1];
            v_fixed = new int[v_num + 1];
            RemovedPoint = new int[v_num + 1];
            AddedPoint = new int[v_num + 1];
            cutPointSet = new int[v_num + 1];
            taburemove = new long[v_num + 1];
            tabuadd = new long[v_num + 1];
            trunkTabuAdd = new long[v_num + 1];
            trunkTabuRemove = new long[v_num + 1];
            child = new int[v_num + 1];
            SF = new int[v_num + 1];
            Stack = new int[v_num + 1];
            onlydominate = new long[v_num + 1];
            onlyDominateByDApos = new long[v_num + 1];

            father = new int[v_num + 1];
            childnum = new int[v_num + 1];
            removedNodeNeighbor = new Array(v_num);
            redundantNodes = new Array(v_num);
            candidateForDApos = new Array(v_num);
            notDominatedByDApos = new Array(v_num);
            constructSolutionHeap = new MyHeap(v_num + 1);
            constructTreeHeap = new MyHeap(v_num + 1);

            /*****************fill_n******************/
            fill_n(pre_deci_step, v_num + 1, v_num);
            fill_n(temp_pre_deci_step, v_num + 1, v_num);
            fill_n(subWeight, v_num + 1, 0);
            fill_n(v_threshold, v_num + 1, 0);
            fill_n(v_degree, v_num + 1, 0);
            fill_n(child, v_num + 1, 0);
            fill_n(taburemove, v_num + 1, 0);
            fill_n(tabuadd, v_num + 1, 0);
            fill_n(trunkTabuAdd, v_num + 1, 0);
            fill_n(trunkTabuRemove, v_num + 1, 0);
            fill_n(v_in_c, v_num + 1, 0);
            fill_n(last_v_in_c, v_num + 1, 0);
            fill_n(lastTurn_v_in_c, v_num + 1, 0);
            fill_n(score, v_num + 1, 0);
            fill_n(subscore, v_num + 1, 0);
            fill_n(score_in_d, v_num + 1, 0);
            fill_n(conf_change, v_num + 1, 1);
            fill_n(time_stamp, v_num + 1, 0);
            fill_n(time_stamp_tree, v_num + 1, 0);
            fill_n(dominated, v_num + 1, 0);
            fill_n(best_dominated, v_num + 1, 0);
            fill_n(v_fixed, v_num + 1, 0);
            fill_n(CV, v_num + 1, 0);
            fill_n(dnf, v_num + 1, 0);
            fill_n(low, v_num + 1, 0);
            fill_n(isCut, v_num + 1, 0);
            fill_n(cutPointSet, v_num + 1, 0);
            fill_n(first, v_num + 1, -1);
            fill_n(dominatedByDApos, v_num + 1, 0);
            fill_n(onlydominate, v_num + 1, 0);
            fill_n(onlyDominateByDApos, v_num + 1, 0);
            /***********************************/
            for (int v = 1; v <= v_num; v++)
                pre[v] = v;
        }
        if (strcmp(tempStr, "v") == 0)
        {
            int v;
            double vertex_weight;
            infile >> v >> vertex_weight;
            frequency[v] = 1;
            weight[v] = vertex_weight;
            totalweight += vertex_weight;
        }
        if (strcmp(tempStr, "e") == 0)
        {
            int u, v;
            infile >> u >> v;
            edge[edgeCount].u = u;
            edge[edgeCount].v = v;
            edge[edgeCount].next = first[edge[edgeCount].u];
            first[edge[edgeCount].u] = edgeCount;

            edge[edgeCount + 1].u = edge[edgeCount].v;
            edge[edgeCount + 1].v = edge[edgeCount].u;
            edge[edgeCount + 1].next = first[edge[edgeCount + 1].u];
            first[edge[edgeCount + 1].u] = edgeCount + 1;

            ++v_degree[u];
            ++v_degree[v];
            edgeCount += 2;
        }
        infile >> tempStr;
    }
    bestWeight = totalweight;
    ave_weight = totalweight / v_num;

    infile.close();
    constructSolutionHeap->registerCmpCallback([](int a, int b) -> bool {
        int scoreA = (int)(ave_weight * subscore[a] / weight[a]);
        int scoreB = (int)(ave_weight * subscore[b] / weight[b]);
        if (isUseBestSolutionInfo && isRecordBestSolution)
        {
            if (best_v_in_c[a] != 1)
                scoreA *= reduceRate;
            if (lastTurn_v_in_c[a] != 1)
                scoreA *= reduceRate;
            if (best_v_in_c[b] != 1)
                scoreB *= reduceRate;
            if (lastTurn_v_in_c[b] != 1)
                scoreB *= reduceRate;
        }
        if (scoreA > scoreB || (scoreA == scoreB && subWeight[a] > subWeight[b]) || (scoreA == scoreB && subWeight[a] == subWeight[b] && a < b))
        {
            return true;
        }
        return false;
    });
    constructSolutionHeap->registerScoreCallback([](int v) -> int {
        return (int)(ave_weight * subscore[v] / weight[v]);
    });
    constructTreeHeap->registerCmpCallback([](int a, int b) -> bool {
        std::function<int(int)> scoreFunc = [](int v) -> int {
            return (int)(ave_weight * score_in_d[v] / weight[v]);
        };
        std::function<int(int)> subscoreFunc = [](int v) -> int {
            return (int)(ave_weight * subscore[v]) / weight[v];
        };
        if (scoreFunc(a) > scoreFunc(b) || (scoreFunc(a) == scoreFunc(b) && subscoreFunc(a) < subscoreFunc(b)) || (scoreFunc(a) == scoreFunc(b) && subscoreFunc(a) == subscoreFunc(b) && subWeight[a] > subWeight[b])){return true;}
        return false;
    });

    v_adj[0] = 0;
    v_edges[0] = 0;
    for (int v = 1; v < v_num + 1; v++)
    {
        v_adj[v] = new long[v_degree[v]];
        v_edges[v] = new int[v_degree[v]];
    }

    int *v_degree_tmp = new int[v_num + 1];
    fill_n(v_degree_tmp, v_num + 1, 0);

    for (int e = 0; e < edgeCount/2; e++)
    {
        int v1 = edge[2 * e].u;
        int v2 = edge[2 * e].v;

        v_edges[v1][v_degree_tmp[v1]] = e;
        v_edges[v2][v_degree_tmp[v2]] = e;

        v_adj[v1][v_degree_tmp[v1]] = v2;
        v_adj[v2][v_degree_tmp[v2]] = v1;

        v_degree_tmp[v1]++;
        v_degree_tmp[v2]++;
    }
    for (int v = 1; v <= v_num; v++)
    {
        // WV[v] = v_degree[v] + 1;
        if (v_degree[v] > v_degree[maxDegreeNode])
            maxDegreeNode = v;
    }
    //    cout<<v_degree[maxDegreeNode];
    averagedegree = e_num / v_num * 2;

    delete[] v_degree_tmp;

    return 0;
}

void FreeMemory()
{
    int v;
    for (v = 0; v < v_num + 1; v++)
    {
        delete[] v_adj[v];
        delete[] v_edges[v];
    }
    delete constructSolutionHeap;
    delete removedNodeNeighbor;
    delete redundantNodes;
    delete greyPointArray;
    delete greyPointArrayDApos;
    delete undomPointArray;
    delete candidateArray;
    delete candidateAposArray;
    delete leafArray;
    delete candidateForDApos;
    delete notDominatedByDApos;
    delete[] dominated;
    delete[] dominatedByDApos;
    delete[] conf_change;
    delete[] best_v_in_c;
    delete[] subWeight;
    delete[] pre_deci_step;
    delete[] temp_pre_deci_step;
    delete[] v_in_c;
    delete[] last_v_in_c;
    delete[] lastTurn_v_in_c;
    delete[] v_degree;
    delete[] v_adj;
    delete[] v_edges;
    delete[] time_stamp;
    delete[] time_stamp_tree;
    delete[] score;
    delete[] pre;
    delete[] CV;
    delete[] WV;
    delete[] S;
    delete[] isInS;
    delete[] RemovedPoint;
    delete[] AddedPoint;
    delete[] cutPointSet;
    delete[] edge;
    delete[] frequency;
    delete[] weight;
    delete[] subscore;
    delete[] score_in_d;
    delete[] taburemove;
    delete[] tabuadd;
    delete[] trunkTabuAdd;
    delete[] trunkTabuRemove;
    delete[] first;
    delete[] child;
    delete[] SF;
    delete[] Stack;
    delete[] v_threshold;
    delete[] onlydominate;
    delete[] fixedSet;
    delete[] father;
    delete[] childnum;
    delete[] isCut;
    delete[] initIndex;
    delete[] index_in_candidate;
    delete[] candidate;
}

#endif
