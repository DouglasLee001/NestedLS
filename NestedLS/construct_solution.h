#ifndef _CONSTRUCT_SOLUTION_H
#define _CONSTRUCT_SOLUTION_H

#include "basic.h"

void resetScore()
{
    fill_n(score, v_num + 1, 0);
    fill_n(subscore, v_num + 1, 0);
    for (size_t i = 1; i < v_num + 1; i++)
    {
        score[i] = v_degree[i] + 1;
        subscore[i] = frequency[i];
        for (size_t j = 0; j < v_degree[i]; j++)
        {
            subscore[i] += frequency[v_adj[i][j]];
        }
    }
}

void resetScoreInD(){
    fill_n(score_in_d, v_num + 1, 0);
    for (size_t v = 0; v < v_num + 1; v++){
        if (v_in_c[v] == 1){
            score_in_d[v]++;
            for (size_t j = 0; j < v_degree[v]; j++){
                int neighbor = v_adj[v][j];
                if (v_in_c[neighbor] == 1){
                    score_in_d[v]++;
                }
            }
        }
    }
}

void constructIncreaseDominate(int v, int v_dominater)
{
    if (dominated[v] == 0)
    {
        --score[v];
        subscore[v] -= frequency[v];
        for (size_t i = 0; i < v_degree[v]; i++)
        {
            int u = v_adj[v][i];
            --score[u];
            subscore[u] -= frequency[v];
            constructSolutionHeap->modify(u);
        }
        Dom(v);
        onlydominate[v] = v_dominater;
        if (v != v_dominater)
        {
            greyPointArray->insert_element(v);
            constructSolutionHeap->insert(v);
        }
    }
    else if (dominated[v] == 1)
    {
        if (v_in_c[v] == 1)
        {
            ++score[v];
            subscore[v] += frequency[v];
        }
        else
        {
            int v_dominater = onlydominate[v];
            ++score[v_dominater];
            subscore[v_dominater] += frequency[v];
        }
    }
    ++dominated[v];
}

void treeIncreaseDominate(int v, int v_dominater, bool isConstruct = true){
    if (dominatedByDApos[v] == 0){
        --score_in_d[v];
        for (size_t i = 0; i < v_degree[v]; i++){
            int u = v_adj[v][i];
            if (v_in_c[u] == 1){
                --score_in_d[u];
                if (isConstruct){constructTreeHeap->modify(u);}
            }
        }
        onlyDominateByDApos[v] = v_dominater;
        notDominatedByDApos->delete_element(v);
        if (isConstruct){
            if (v != v_dominater){
                constructTreeHeap->insert(v);
                greyPointArrayDApos->insert_element(v);
            }
        }
        else
            greyPointArrayDApos->insert_element(v);
    }
    else if (dominatedByDApos[v] == 1){
        if (candidateAposArray->is_in_array(v)){++score_in_d[v];}
        else{
            int v_dominater = onlyDominateByDApos[v];
            ++score_in_d[v_dominater];
            }
    }
    ++dominatedByDApos[v];
}

void treeDecreaseDominate(int v, bool isConstruct = false){
    if (dominatedByDApos[v] == 1){
        ++score_in_d[v];
        for (size_t i = 0; i < v_degree[v]; i++){
            int neighbor = v_adj[v][i];
            if (v_in_c[neighbor] == 1){
                ++score_in_d[neighbor];
                if (isConstruct){constructTreeHeap->modify(neighbor);}
            }
        }
        greyPointArrayDApos->delete_element(v);
        notDominatedByDApos->insert_element(v);
    }
    else if (dominatedByDApos[v] == 2){
        if (candidateAposArray->is_in_array(v)){--score_in_d[v];}
        else{
            for (size_t i = 0; i < v_degree[v]; i++){
                long u = v_adj[v][i];
                if (candidateAposArray->is_in_array(u)){
                    --score_in_d[u];
                    onlyDominateByDApos[v] = u;
                    break;
                }
            }
        }
    }
    --dominatedByDApos[v];
}

void constructAdd(int v, bool isInitial = false)
{
    if (v_in_c[v] == 1)
        return;
    v_in_c[v] = 1;
    if (isInitial){last_v_in_c[v] = 1;}
    c_size++;
    currentWeight += weight[v];
    greyPointArray->delete_element(v);

    if (currentMode == ChooseMode::ModeC)
    {
        pre_deci_step[v] = add_step++;
    }

    int new_score = -score[v];
    int new_subscore = -subscore[v];
    constructIncreaseDominate(v, v);
    for (size_t i = 0; i < v_degree[v]; i++)
    {
        int u = v_adj[v][i];
        constructIncreaseDominate(u, v);
    }
    score[v] = new_score;
    subscore[v] = new_subscore;
}

void prepareForLocalSearch()
{
    ConstructSolution();
    if(TimeElapsed()>cutoff_time)return;
    bestWeight = currentWeight;
    candidateArray->clear();
    for (int i = 1; i <= v_num; ++i){candidateArray->insert_element(i);}
    MarkCut();
    for (int i = 0; i < cutIndex; ++i){
        if (isCut[cutPointSet[i]] != 0){
            int cutPoint = cutPointSet[i];
            if (v_fixed[cutPoint] == 0){
                v_fixed[cutPoint] = 1;
                fixedSet[fixedNum++] = cutPoint;
            }
        }
    }
    ResetCandidate();
    RemoveRedundant();
    UpdateBestSolution();
}

int NewSolutionChooseVFromMethodA()
{
    int best_add_v = -1;
    int cscore;
    int best_cscore = -weightthreshold;
    if (c_size == 0){
        for (size_t i = 1; i < v_num + 1; i++){
            cscore = (int)(ave_weight * subscore[i] / weight[i]);
            if (cscore > best_cscore){
                best_add_v = i;
                best_cscore = cscore;
            }
            else if (cscore == best_cscore)
            {
                if (subWeight[i] > subWeight[best_add_v])
                {
                    best_add_v = i;
                }
            }
        }
    }
    else
    {
        best_add_v = constructSolutionHeap->remove_first();
    }
    return best_add_v;
}
int NewSolutionChooseVFromMethodB()
{
    return -1;
}
int NewSolutionChooseVFromMethodC()
{
    int best_add_v = -1;
    int cscore;
    int best_cscore = -weightthreshold;
    if (c_size == 0)
    {
        for (size_t i = 1; i < v_num + 1; i++)
        {
            cscore = (int)(ave_weight * subscore[i] / weight[i]);
            if (cscore > best_cscore)
            {
                best_add_v = i;
                best_cscore = cscore;
            }
            else if (cscore == best_cscore)
            {
                if (pre_deci_step[i] > pre_deci_step[best_add_v])
                {
                    best_add_v = i;
                }
            }
        }
    }
    else
    {
        for (size_t i = 0; i < greyPointArray->size(); i++)
        {
            int greyPoint = greyPointArray->element_at(i);
            cscore = (int)(ave_weight * subscore[greyPoint] / weight[greyPoint]);
            if (cscore > best_cscore)
            {
                best_add_v = greyPoint;
                best_cscore = cscore;
            }
            else if (cscore == best_cscore)
            {
                if (pre_deci_step[greyPoint] > pre_deci_step[best_add_v])
                {
                    best_add_v = greyPoint;
                }
            }
        }
    }
    return best_add_v;
}
int NewSolutionChooseVFromMethodD()
{
    return -1;
}

void ConstructSolution()
{
    undomPointArray->clear();
    for (int v = 1; v < v_num + 1; ++v)
    {
        undomPointArray->insert_element(v);
    }
    constructSolutionHeap->clear();
    resetScore();
    int construct_step=0;
    while (undomPointArray->size() != 0)
    {
        if(construct_step++>1000){
            construct_step=0;
            if(TimeElapsed()>cutoff_time)return;
        }
        int addV = -1;
        switch (currentMode)
        {
        case ChooseMode::ModeA:
            addV = NewSolutionChooseVFromMethodA();
            break;
        case ChooseMode::ModeB:
            addV = NewSolutionChooseVFromMethodA();
            break;
        case ChooseMode::ModeC:
            addV = NewSolutionChooseVFromMethodC();
            break;
        case ChooseMode::ModeD:
            addV = NewSolutionChooseVFromMethodD();
            break;
        }
        if (addV != -1)
        {
            constructAdd(addV);
            if (currentMode == ChooseMode::ModeC)
            {
                temp_pre_deci_step[addV] = add_step++;
            }
        }
    }
}

void ConstructSolutionToBestWeight()
{
    undomPointArray->clear();
    for (int v = 1; v < v_num + 1; ++v)
    {
        undomPointArray->insert_element(v);
    }
    constructSolutionHeap->clear();
    resetScore();

    while (undomPointArray->size() != 0 && currentWeight < bestWeight)
    {
        int addV = -1;
        switch (currentMode)
        {
        case ChooseMode::ModeA:
            addV = NewSolutionChooseVFromMethodA();
            break;
        case ChooseMode::ModeB:
            addV = NewSolutionChooseVFromMethodA();
            break;
        case ChooseMode::ModeC:
            addV = NewSolutionChooseVFromMethodC();
            break;
        case ChooseMode::ModeD:
            addV = NewSolutionChooseVFromMethodD();
            break;
        }
        if (addV != -1)
        {
            if (currentWeight + weight[addV] > bestWeight)
            {
                return;
            }
            constructAdd(addV);
            if (currentMode == ChooseMode::ModeC)
            {
                temp_pre_deci_step[addV] = add_step++;
            }
        }
    }
    minUndom = undomPointArray->size();
}

double getNewSolutionRepeatDegree()
{
    int lastCount = 0;
    int newCount = 0;
    int repeatCount = 0;
    for (size_t v = 1; v < v_num + 1; v++)
    {
        if (v_in_c[v] == 1 && last_v_in_c[v] == 1)
        {
            repeatCount++;
        }
        else if (v_in_c[v] == 1)
        {
            newCount++;
        }
        else if (last_v_in_c[v] == 1)
        {
            lastCount++;
        }
    }
    return repeatCount * 1.0 / (lastCount + repeatCount + newCount) * 100;
}

#endif
