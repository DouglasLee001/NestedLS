
#include "basic.h"
#include <queue>

void ResetCandidate()
{
    candidateArray->clear();
    for (int v = 1; v < v_num + 1; v++)
    {
        if (v_in_c[v] == 1 && v_fixed[v] == 0)
        {
            candidateArray->insert_element(v);
        }
        else
        {
            if (dominated[v] > 0)
            {
                greyPointArray->insert_element(v);
            }
        }
    }
}

void updateSolution()
{
    if (undomPointArray->size() <= minUndom)
    {
        minUndom = undomPointArray->size();
        if (isUseBestSolutionInfo)
        {
            memcpy(lastTurn_v_in_c, v_in_c, sizeof(int) * (v_num + 1));
        }
        if (minUndom == 0)
        {
            best_c_size = c_size;
            best_comp_time = TimeElapsed();
            best_step = step;
            bestWeight = currentWeight;
            if (isUseBestSolutionInfo)
            {
                isRecordBestSolution = true;
                memcpy(best_v_in_c, v_in_c, sizeof(int) * (v_num + 1));
            }
//            cout << "weight: " << currentWeight << "\ttime: " << best_comp_time << "\tstepAction:" << stepAction <<"\tisRecord:"<<(isRecordBestSolution ? "true":"false ")<<endl;
        }
    }
}

void UpdateBestSolution()
{
    if ((int)currentWeight < (int)bestWeight)
    {
        best_c_size = c_size;
        best_comp_time = TimeElapsed();
        best_step = step;
        bestWeight = currentWeight;
        if (isUseBestSolutionInfo)
        {
            isRecordBestSolution = true;
            memcpy(best_v_in_c, v_in_c, sizeof(int) * (v_num + 1));
        }

//        cout << "weight: " << currentWeight << "\ttime: " << best_comp_time << "\tstepAction:" << stepAction <<"\tisRecord:"<<(isRecordBestSolution ? "true":"false ")<<endl;
    }
}

void cutPoint(int cur, int fa)
{
    int child = 0;
    dnf[cur] = low[cur] = ++ind;
    for (int i = 0; i < v_degree[cur]; i++)
    {
        int u = v_adj[cur][i];
        if (v_in_c[u] == 1)
        {
            if (dnf[u] == 0)
            {
                child++;
                cutPoint(u, cur);
                low[cur] = min(low[cur], low[u]);
                if ((cur != root && low[u] >= dnf[cur]) || (cur == root && child >= 2))
                {
                    isCut[cur] = 1;
                    cutPointSet[cutIndex++] = cur;
                }
            }
            else if (u != fa)
            {
                low[cur] = min(low[cur], dnf[u]);
            }
        }
    }
    if (score[cur] > maxScore)
        maxPoint = cur;
}

void cutPointNoRecur(int root)
{
    int k, u, v, inde, top = 0;
    inde = 0;
    for (u = 0; u < candidateArray->size(); u++)
    {
        v = candidateArray->element_at(u);
        SF[v] = first[v];
    }
    for (u = 0; u < fixedNum; u++)
        SF[fixedSet[u]] = first[fixedSet[u]];
    Stack[top] = root;
    dnf[root] = low[root] = ++inde;
    while (top >= 0)
    {
        k = SF[Stack[top]];

        if (k != -1)
        {
            if (v_in_c[edge[k].v] == 0)
            {
                SF[Stack[top]] = edge[k].next;
                continue;
            }
            SF[Stack[top]] = edge[k].next;
            if (dnf[edge[k].v] == 0)
            {
                child[Stack[top]]++;
                Stack[++top] = edge[k].v;
                low[edge[k].v] = dnf[edge[k].v] = ++inde;
            }
            else
            {
                low[Stack[top]] = (low[Stack[top]] < dnf[edge[k].v]) ? low[Stack[top]] : dnf[edge[k].v];
            }
        }
        else
        {
            if (top > 0)
            {
                u = Stack[top - 1];
                v = Stack[top];
                low[u] = (low[u] < low[v]) ? low[u] : low[v];
                if ((u != root && low[v] >= dnf[u]) || (u == root && child[u] >= 2))
                {
                    isCut[u] = 1;
                    cutPointSet[cutIndex++] = u;
                }
            }
            top--;
        }
    }
}

void cutPoint1(long root)
{
    queue<int> qu;
    qu.push(root);
    dnf[root] = ++ind;
    while (!qu.empty())
    {
        int cur = qu.front();
        childnum[cur] = 0;
        qu.pop();
        int cur_degree=v_degree[cur];
        for (int i = 0; i < v_degree[cur]; i++)
        {
            int u = v_adj[cur][i];
            if (v_in_c[u] == 1 && dnf[u] == 0)
            {
                dnf[u] = ++ind;
                childnum[cur]++;
                qu.push(u);
                father[u] = cur;
            }
        }
        if ((v_fixed[cur] == 0) && (childnum[cur] == 0 || (childnum[cur] == 1 && cur == root)))
        {
            leafArray->insert_element(cur);
        }
    }
}

bool judgeCut(int node)
{
    root = node;
    cutPoint(node, node);
    return (isCut[node] == 1);
}
void updateWeight(){
    printDebugMsg("updateWeight");
    double avg_weight = totalweight / v_num;
    for (int i = 1; i <= v_num; i++){
        int newweight = weightAdjustCoefficient * frequency[i] + avg_weight * (1 - weightAdjustCoefficient);
        if (newweight < 1){newweight = 1;}
        int tobeminus = frequency[i] - newweight;
        minusWeight(i, tobeminus);
        totalweight -= tobeminus;
    }
}

void updateSubscore()
{
    fill_n(frequency, v_num + 1, 1);
    memcpy(subscore, score, sizeof(int) * (v_num + 1));
    totalweight = v_num;
}
