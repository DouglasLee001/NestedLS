#ifndef _FASTDS_PURE_H
#define _FASTDS_PURE_H

#include "basic.h"

inline void Undom(int v){undomPointArray->insert_element(v);}


inline void Dom(int v){undomPointArray->delete_element(v);}

void updateRedundantV(int v){
    if (!moduleRemoveRedundant){return;}
    if (v_in_c[v] == 1 && score[v] == 0 && v_fixed[v] == 0) {redundantNodes->insert_element(v);}
    else{redundantNodes->delete_element(v);}
}

void RemoveRedundant(){
    if (!moduleRemoveRedundant){return;}
    MarkCut();
    while (redundantNodes->size() != 0){
        int bestRemoveV = -1;
        double bestRemoveVWeight = 0.0;
        for (size_t i = 0; i < redundantNodes->size(); i++){
            int redundantV = redundantNodes->element_at(i);
            if (weight[redundantV] > bestRemoveVWeight && isCut[redundantV] == 0 && v_in_c[redundantV] == 1 && v_fixed[redundantV] == 0){
                bestRemoveV = redundantV;
                bestRemoveVWeight = weight[redundantV];
            }
        }
        if (bestRemoveV != -1){
            Remove(bestRemoveV);
            MarkCut();
        }
        else{break;}
    }
}

void Restart(){
    stepAction = 1;
    NOimprovementstep = 0;
    Noimprovementstepfocus = 0;
    greyPointArray->clear();
    leafArray->clear();
    fill_n(v_in_c, v_num + 1, 0);
    fill_n(dominated, v_num + 1, 0);
    step += 20;
    fill_n(conf_change, v_num + 1, 1);
    fill_n(isCut, v_num + 1, 0);
    fill_n(cutPointSet, v_num + 1, 0);
    c_size = 0;
    cutIndex = 0;
    currentWeight = 0;
    if (currentMode == ChooseMode::ModeC){add_step = 0;}
    if (lastRepeatDegree > restartUpdateFrequencyThreshold){updateWeight();}
    ConstructSolutionToBestWeight();
    if (currentMode == ChooseMode::ModeA || currentMode == ChooseMode::ModeB){fill_n(subWeight, v_num + 1, 0);}
    else if (currentMode == ChooseMode::ModeC){
        memcpy(pre_deci_step, temp_pre_deci_step, sizeof(int) * (v_num + 1));
        fill_n(temp_pre_deci_step, v_num + 1, v_num);
    }
    ResetCandidate();
    lastRepeatDegree = getNewSolutionRepeatDegree();
    memcpy(last_v_in_c, v_in_c, sizeof(int) * (v_num + 1));
    if (undomPointArray->empty() && currentWeight < bestWeight){UpdateBestSolution();}
    MarkCut();
}

void increase_dominate(long v, long source_v){
    if (dominated[v] == 0){
        --score[v];
        subscore[v] -= frequency[v];
        updateRedundantV(v);
        for (int i = 0; i < v_degree[v]; ++i){
            long u = v_adj[v][i];
            --score[u];
            subscore[u] -= frequency[v];
            updateRedundantV(v);
        }
        Dom(v);
        onlydominate[v] = source_v;
        greyPointArray->insert_element(v);
    }
    else if (dominated[v] == 1){
        if (v_in_c[v] == 1){
            ++score[v];
            subscore[v] += frequency[v];
            updateRedundantV(v);
        }
        else{long u = onlydominate[v];
            ++score[u];
            subscore[u] += frequency[v];
            updateRedundantV(u);
        }
    }
    ++dominated[v];
}

void addUpdate(int v, bool istrunk = false){
    int u = -1;
    int uinc;
    childnum[v] = 0;
    for (int i = 0; i < v_degree[v]; i++){
        u = v_adj[v][i];
        if (!istrunk){
            if (v_in_c[u] == 1){
                uinc = u;
                if (childnum[u] != 0 && u != root){
                    childnum[u]++;
                    father[v] = u;
                    if (v_fixed[v] == 0)
                        leafArray->insert_element(v);
                    return;
                }
            }
        }
        else{
            if (candidateAposArray->is_in_array(u)){
                uinc = u;
                if (childnum[u] != 0 && u != root){
                    childnum[u]++;
                    father[v] = u;
                    if (v_fixed[v] == 0){leafArray->insert_element(v);}
                    return;
                }
            }
        }
    }
    if (uinc == root){
        if (leafArray->is_in_array(uinc)){
            if (v_fixed[v] == 0){leafArray->insert_element(v);}
            leafArray->delete_element(uinc);
        }
        else{
            if (v_fixed[v] == 0) {leafArray->insert_element(v);}
        }
        father[uinc] = v;
        childnum[v] = 1;
        father[v] = 0;
        root = v;
    }
    else{
        childnum[uinc]++;
        father[v] = uinc;
        if (leafArray->is_in_array(uinc)){
            if (v_fixed[v] == 0){leafArray->insert_element(v);}
            leafArray->delete_element(uinc);
        }
        else{
            if (v_fixed[v] == 0){leafArray->insert_element(v);}
        }
    }
}


bool Add(int v){
    if (v_in_c[v] == 1 || v < 0){return false;}
    greyPointArray->delete_element(v);
    if (currentMode == ChooseMode::ModeA){subWeight[v]++;}
    int new_score = -score[v];
    int new_subscore = -subscore[v];
    increase_dominate(v, v);
    for (int i = 0; i < v_degree[v]; ++i){
        long u = v_adj[v][i];
        increase_dominate(u, v);
    }
    score[v] = new_score;
    subscore[v] = new_subscore;
    v_in_c[v] = 1;
    updateRedundantV(v);
    currentWeight += weight[v];
    c_size++;
    candidateArray->insert_element(v);
    int u;
    int tmp = conf_change[v];
    for (int i = 0; i < v_degree[v]; i++)
    {
        u = v_adj[v][i];
        conf_change[u] = 1; //
    }
    conf_change[v] = tmp;
    taburemove[v] = (step + 1 + rand() % 3);
    return true;
}


void addToDApos(int v, bool isConstruct = true){
    if (candidateAposArray->is_in_array(v)){return;}
    greyPointArrayDApos->delete_element(v);
    int new_score_in_d = -score_in_d[v];
    treeIncreaseDominate(v, v, isConstruct);
    for (size_t i = 0; i < v_degree[v]; i++){
        int neighbor = v_adj[v][i];
        if (v_in_c[neighbor] == 1){treeIncreaseDominate(neighbor, v, isConstruct);}
    }
    score_in_d[v] = new_score_in_d;
    candidateAposArray->insert_element(v);
    candidateForDApos->delete_element(v);
    curTrunkWeight += weight[v];
    if (!isConstruct){trunkTabuRemove[v] = trunkStep + (1 + rand() % 3);}
}

void removeFromDApos(int v, bool isConstruct = false)
{
    if (!candidateAposArray->is_in_array(v)){return;}
    greyPointArrayDApos->insert_element(v);
    int new_score_in_d = -score_in_d[v];
    treeDecreaseDominate(v, isConstruct);
    for (size_t i = 0; i < v_degree[v]; i++){
        int neighbor = v_adj[v][i];
        if (v_in_c[neighbor] == 1){treeDecreaseDominate(neighbor, isConstruct);}
    }
    score_in_d[v] = new_score_in_d;
    candidateAposArray->delete_element(v);
    candidateForDApos->insert_element(v);
    curTrunkWeight -= weight[v];
    if (!isConstruct){trunkTabuAdd[v] = trunkStep + (1 + rand() % 3);}
}

void decrease_dominate(int v){
    if (dominated[v] == 1){
        ++score[v];
        subscore[v] += frequency[v];
        updateRedundantV(v);
        for (int i = 0; i < v_degree[v]; ++i)
        {
            int u = v_adj[v][i];
            ++score[u];
            subscore[u] += frequency[v];
            updateRedundantV(u);
        }
        Undom(v);
        greyPointArray->delete_element(v);
    }
    else if (dominated[v] == 2)
    {
        if (v_in_c[v]){
            --score[v];
            subscore[v] -= frequency[v];
            updateRedundantV(v);
        }
        else{
            for (int i = 0; i < v_degree[v]; ++i){
                long u = v_adj[v][i];
                if (v_in_c[u]){
                    --score[u];
                    subscore[u] -= frequency[v];
                    onlydominate[v] = u;
                    break;
                }
            }
        }
    }
    --dominated[v];
}

void removeUpdate(int v, bool istrunk = false){
    if (v != root){
        int fa = father[v];
        childnum[fa]--;
        if (fa != root){
            if (childnum[fa] == 0 && v_fixed[fa] == 0){
                leafArray->insert_element(fa);
                leafArray->delete_element(v);
            }
            else{leafArray->delete_element(v);}
        }
        else{
            if (childnum[fa] == 1 && v_fixed[fa] == 0){
                leafArray->insert_element(fa);
                leafArray->delete_element(v);
            }
            else{leafArray->delete_element(v);}
        }
    }
    else{
        int chi;
        for (int i = 0; i < v_degree[v]; i++){
            int u = v_adj[v][i];
            if (istrunk == false){
                if (v_in_c[u] == 1 && father[u] == v){
                    chi = u;
                    break;
                }
            }
            else{
                if (candidateAposArray->is_in_array(u) == 1 && father[u] == v){
                    chi = u;
                    break;
                }
            }
        }
        root = chi;
        father[root] = 0;
        if (childnum[chi] == 1 && v_fixed[chi] == 0){
            leafArray->insert_element(chi);
            leafArray->delete_element(v);
        }
        else{leafArray->delete_element(v);}
    }
    father[v] = 0;
    childnum[v] = 0;
}

bool Remove(int v){
    greyPointArray->insert_element(v);
    v_in_c[v] = 0;
    currentWeight -= weight[v];
    c_size--;
    candidateArray->delete_element(v);
    int new_score = -score[v];
    int new_subscore = -subscore[v];
    decrease_dominate(v);
    int neighbours = v_degree[v];
    for (int i = 0; i < neighbours; ++i){
        int u = v_adj[v][i];
        decrease_dominate(u);
    }
    score[v] = new_score;
    subscore[v] = new_subscore;
    updateRedundantV(v);
    int u;
    for (int i = 0; i < v_degree[v]; i++){
        u = v_adj[v][i];
        conf_change[u] = 1;
    }
    conf_change[v] = 0;
    tabuadd[v] = step + 1;
    return true;
}


void addWeight(int node){
    int increment = 5;
    frequency[node] += increment;
    for (int i = 0; i < v_degree[node]; i++)
        subscore[v_adj[node][i]] += increment;
    subscore[node] += increment;
    if (currentMode == ChooseMode::ModeB){subWeight[node]++;}
}

void minusWeight(int node, int tobeminus){
    frequency[node] -= tobeminus;
    if (dominated[node] == 0){
        for (int i = 0; i < v_degree[node]; i++){subscore[v_adj[node][i]] -= tobeminus;}
        subscore[node] -= tobeminus;
    }
    else if (v_in_c[node] == 0){
        if (dominated[node] == 1){
            long u = onlydominate[node];
            subscore[u] += tobeminus;
        }
    }
}

bool cmp(int a, int b){
    return (subscore[a] / weight[a] > subscore[b] / weight[b]);
}



int ChooseRemoveVTopofBMS(const int count, int choice){
    if (candidateArray->empty()){return -1;}
    int v, i;
    int best_score = -weightthreshold;
    int best_remove_v = -1;
    vector<int> toberemoved1(count);
    int topIndex = 0;
    int cscore;
    for (i = 0; i < count; i++){
        if (choice == 1){
            v = leafArray->rand_element();
            if (v_fixed[v] == 1)
                continue;
        }
        else{
            v = candidateArray->rand_element();
            if (v_fixed[v] == 1 || isCut[v] == 1)
                continue;
        }
        cscore = (int)(ave_weight * subscore[v] / weight[v]);
        toberemoved1[topIndex++] = v;
        if (step > taburemove[v]){
            if (cscore > best_score){
                best_remove_v = v;
                best_score = cscore;
            }
            else if (cscore == best_score){
                if (((int)(ave_weight * dominated[v] / weight[v]) < (int)(ave_weight * dominated[best_remove_v] / weight[best_remove_v])) || ((int)(ave_weight * dominated[v] / weight[v]) == (int)(ave_weight * dominated[best_remove_v] / weight[best_remove_v]) && time_stamp[v] < time_stamp[best_remove_v]))
                    best_remove_v = v;
            }
        }
    }
    if (best_remove_v != -1){return best_remove_v;}
    else{
        if (topIndex > 0){return toberemoved1[rand() % topIndex];}
        else{return -1;}
    }
}

int ChooseRemoveVFromArray(Array *removedNodeNeighbor, int choice){
    int best_remove_v = -1;
    int cscore;
    int best_cscore = -weightthreshold;
    for (size_t i = 0; i < removedNodeNeighbor->size(); i++){
        int v = removedNodeNeighbor->element_at(i);
        if (choice == 0){
            if (v_in_c[v] == 0 || v_fixed[v] == 1 || isCut[v] == 1) {continue;}
        }
        else{
            if (v_in_c[v] == 0 || !leafArray->is_in_array(v) || v_fixed[v] == 1){continue;}
        }
        cscore = (int)(ave_weight * subscore[v] / weight[v]);
        if (step > taburemove[v]){
            if (cscore > best_cscore){
                best_remove_v = v;
                best_cscore = cscore;
            }
            else if (cscore == best_cscore){
                if (((int)(ave_weight * dominated[v] / weight[v]) < (int)(ave_weight * dominated[best_remove_v] / weight[best_remove_v])) || ((int)(ave_weight * dominated[v] / weight[v]) == (int)(ave_weight * dominated[best_remove_v] / weight[best_remove_v]) && time_stamp[v] < time_stamp[best_remove_v])){best_remove_v = v;}
            }
        }
    }
    return best_remove_v;
}

int ChooseRemoveVFromDApos(const int count = 50){
    int best_remove_v = -1;
    int best_score = -weightthreshold;
    int cscore;
    vector<int> toberemoved(count);
    int index = 0;
    int v;
    if (!candidateAposArray->empty()){
        for (size_t i = 0; i < count; i++){
            v = leafArray->rand_element();
            if (v_fixed[v] == 1){continue;}
            cscore = (int)(ave_weight * score_in_d[v] / weight[v]);
            toberemoved[index++] = v;
            if (trunkStep > trunkTabuRemove[v]){
                if (cscore > best_score){
                    best_remove_v = v;
                    best_score = cscore;
                }
                else if (cscore == best_score){
                    int subscore1_best = (int)(ave_weight * subscore[best_remove_v]) / weight[best_remove_v];
                    int subscore1_v = (int)(ave_weight * subscore[v]) / weight[v];
                    if (subscore1_v > subscore1_best){best_remove_v = v;}
                    else if (subscore1_best == subscore1_v){
                        if (subWeight[v] < subWeight[best_remove_v]){best_remove_v = v;}
                    }
                }
            }
        }
    }
    if (best_remove_v != -1){return best_remove_v;}
    else{
        if (index > 0){return toberemoved[rand() % index];}
    }
    return best_remove_v;
}


int ChooseAddVsubscorefast(int cnt, int choice = 1){
    int base_v, add_v;
    int cscore;
    int best_score = -weightthreshold;
    int best_add_v = -1;
    map<int, int> m;
    int size = (undomPointArray->size() < cnt) ? undomPointArray->size() : cnt;
    bool useBMS = (undomPointArray->size() <= cnt) ? false : true;
    vector<int> tobeadd1(100);
    int topIndex = 0;
    for (int i = 0; i < size; ++i){
        if (topIndex >= 99){break;}
        base_v = useBMS ? undomPointArray->element_at(rand() % undomPointArray->size()) : undomPointArray->element_at(i);
        for (int j = 0; j < v_degree[base_v]; ++j){
            if (topIndex >= 99){break;}
            add_v = v_adj[base_v][j];
            if (m.find(add_v) == m.end() && greyPointArray->is_in_array(add_v) && (weight[add_v] + currentWeight < bestWeight)){
                m[add_v] = topIndex;
                tobeadd1[topIndex++] = add_v;
                cscore = (int)(ave_weight * subscore[add_v] / weight[add_v]);
                if ((choice == 0 && conf_change[add_v]) == 1 || (choice == 1 && step >= tabuadd[add_v])){
                    if (cscore > best_score){
                        best_add_v = add_v;
                        best_score = cscore;
                    }
                    else if (cscore == best_score){
                        if (((int)(ave_weight * dominated[add_v] / weight[add_v]) > (int)(ave_weight * dominated[best_add_v] / weight[best_add_v])) || ((int)(ave_weight * dominated[add_v] / weight[add_v]) == (int)(ave_weight * dominated[best_add_v] / weight[best_add_v]) && time_stamp[add_v] < time_stamp[best_add_v])){best_add_v = add_v;}
                    }
                }
            }
        }
    }
    if (best_add_v != -1){return best_add_v;}
    else{
        if (topIndex != 0){
            sort(tobeadd1.begin(), tobeadd1.begin() + topIndex, cmp);
            int topof = (topIndex * 0.4) + 1;
            return tobeadd1[rand() % topof];
        }
        else{return -1;}
    }
}

int ChooseAddVForConstructTrunk(const int cnt = 50)
{
    int base_v, add_v;
    int cscore;
    int best_score = -weightthreshold;
    int best_add_v = -1;
    vector<int> tobeadd(cnt);
    map<int, int> m;
    int index = 0;
    for (size_t i = 0; i < cnt; i++){
        if (index >= cnt - 1) break;
        base_v = notDominatedByDApos->rand_element();
        for (size_t j = 0; j < v_degree[base_v]; j++){
            if (index >= cnt - 1) break;
            add_v = v_adj[base_v][j];
            if (v_in_c[add_v] != 1){continue;}
            if (greyPointArrayDApos->is_in_array(add_v)){
                if (m.find(add_v) == m.end()){
                    m[add_v] = index;
                    tobeadd[index++] = add_v;
                    cscore = (int)(ave_weight * score_in_d[add_v] / weight[add_v]);
                    if (trunkStep > trunkTabuAdd[add_v]){
                        if (cscore > best_score){
                            best_add_v = add_v;
                            best_score = cscore;
                        }
                        else if (cscore == best_score){
                            int subscore1_best = (int)(ave_weight * subscore[best_add_v]) / weight[best_add_v];
                            int subscore1_v = (int)(ave_weight * subscore[add_v]) / weight[add_v];
                            if (subscore1_v < subscore1_best){best_add_v = add_v;}
                            else if (subscore1_best == subscore1_v){
                                if (subWeight[add_v] > subWeight[best_add_v]){best_add_v = add_v;}
                            }
                        }
                    }
                }
            }
        }
    }
    if (best_add_v != -1){return best_add_v;}
    else{
        if (index > 0){return tobeadd[rand() % index];}
    }
    return best_add_v;
}

void MarkCut(){
    ind = 0;
    for (int i = 0; i < cutIndex; i++)
        isCut[cutPointSet[i]] = 0;
    leafArray->clear();
    cutIndex = 0;
    if (candidateArray->empty()) {return;}
    root = candidateArray->rand_element();
    cutPointNoRecur(root);
    for (int i = 0; i < candidateArray->size(); i++){
        int v = candidateArray->element_at(i);
        child[v] = low[v] = dnf[v] = 0;
    }
    for (int i = 0; i < fixedNum; i++) {child[fixedSet[i]] = low[fixedSet[i]] = dnf[fixedSet[i]] = 0;}
}


void MarkCuttree(){
    ind = 0;
    leafArray->clear();
    fill_n(father, v_num + 1, 0);
    fill_n(childnum, v_num + 1, 0);
    if (candidateArray->size() > 0){
        root = candidateArray->rand_element();
        cutPoint1(root);
        for (int i = 0; i < candidateArray->size(); i++)
            dnf[candidateArray->element_at(i)] = 0;
        for (int i = 0; i < fixedNum; i++)
            dnf[fixedSet[i]] = 0;
    }
}

void MarkCutTreeTrunk(){
    if (candidateArray->size() > 0){
        fill_n(dominatedByDApos, v_num + 1, 0);
        fill_n(father, v_num + 1, 0);
        fill_n(childnum, v_num + 1, 0);
        fill_n(onlyDominateByDApos, v_num + 1, 0);
        leafArray->clear();
        constructTreeHeap->clear();
        candidateAposArray->clear();
        candidateForDApos->clear();
        notDominatedByDApos->clear();
        greyPointArrayDApos->clear();
        for (int v = 1; v <= v_num; v++){
            if (v_in_c[v] == 1){
                candidateForDApos->insert_element(v);
                notDominatedByDApos->insert_element(v);
            }
        }
        resetScoreInD();
        curTrunkWeight = 0;
        bestTrunkWeight = 0;
        root = candidateForDApos->rand_element();
        constructTreeHeap->insert(root);
        notDominatedByDApos->delete_element(root);
        while (!notDominatedByDApos->empty()){
            int add_to_trunk_v = constructTreeHeap->remove_first();
            if (add_to_trunk_v != root) {addUpdate(add_to_trunk_v, true);}
            addToDApos(add_to_trunk_v);
        }
        bestTrunkWeight = curTrunkWeight;
         localSearchTrunk();
        for (size_t i = 0; i < candidateForDApos->size(); i++){
            int v = candidateForDApos->element_at(i);
            addUpdate(v, true);
        }
        candidateForDApos->clear();
    }
}

void localSearchTrunk(){
    trunkStep += 5;
    try_step = 1000;
    trunkStepAction = 1;
    long TreeNoimprovementstep = 0;
    long improvementCount = 0;
    int neighborSize = 1;
    while (true){
        if (leafArray->empty()){return;}
        if (trunkStepAction % try_step == 0){
            int timenow = TimeElapsed();
            if (timenow > cutoff_time){return;}
        }
        if (trunkStepAction > maxTrunkStepAction || TreeNoimprovementstep > noImproveTrunkTryStep){return;}
        for (size_t i = 0; i < neighborSize; i++){
            int best_removed_v = -1;
            if (leafArray->empty()){break;}
            best_removed_v = ChooseRemoveVFromDApos();
            if (best_removed_v != -1){
                removeFromDApos(best_removed_v);
                removeUpdate(best_removed_v, true);
            }
        }
        while (!notDominatedByDApos->empty()){
            int best_add_v = ChooseAddVForConstructTrunk();
            if (best_add_v != -1){
                addToDApos(best_add_v, false);
                addUpdate(best_add_v, true);
            }
            else{break;}
        }
        if (curTrunkWeight < bestTrunkWeight){
            bestTrunkWeight = curTrunkWeight;
            neighborSize = 1;
            improvementCount++;
        }
        else{
            TreeNoimprovementstep++;
            neighborSize = (neighborSize > maxNeighborSize || neighborSize >= candidateAposArray->size()) ? 1 : neighborSize + 1;
        }
        trunkStepAction++;
        trunkStep++;
    }
}

void localSearchFramework2(){
    while (TimeElapsed() < cutoff_time){
        Framework2CutTree();
        if (isShutDownSearchStep){
            isShutDownSearchStep = false;
            Restart();
            tarjanLoopCount = 1;
            continue;
        }
        Framework2TarjanFocus();
    }
}

void modifyTarjanLoopCount(int startWeight, int stepAction){
    if ((startWeight - currentWeight) * 10000 / stepAction > improveRate){
        if (tarjanLoopCount > tarjanCntLB)
            tarjanLoopCount /= 2;
    }
    else if (tarjanLoopCount <= tarjanCntUB){
        isUseBestSolutionInfo = true;
        tarjanLoopCount *= 2;
    }
    if (tarjanLoopCount > tarjanCntUB)
        isShutDownSearchStep = true;
}

void chooseMarkCuttree(){
    if(c_size>500000) MarkCuttree();
    else MarkCutTreeTrunk();
}

void Framework2CutTree()
{
    int startWeight = currentWeight;
    stepAction = 1;
    try_step = 10000;
    long NoImpro_trystep = noImproveTryStep;
    long TreeNoimprovementstep = 0;
    long improvementCount = 0;
    int neighborSize = 1;
    chooseMarkCuttree();
    while (true){
        if (undomPointArray->empty() && leafArray->empty()){return;}
        if (stepAction % try_step == 0){
            int timenow = TimeElapsed();
            if (timenow > cutoff_time){return;}
            if (stepAction >= gap1){
                if (improvementCount > 10){
                    instance1 -= floor1;
                    if (instance1 < floor1){instance1 = floor1;}
                }
                modifyTarjanLoopCount(startWeight, stepAction);
                return;
            }
            if (TreeNoimprovementstep >= instance1){
                instance1 += floor1;
                if (instance1 > gap1){instance1 = gap1;}
                modifyTarjanLoopCount(startWeight, stepAction);
                return;
            }
        }
        for (size_t i = 0; i < neighborSize; i++){
            int best_removed_v = -1;
            if (leafArray->empty()){
                MarkCut();
                best_removed_v = ChooseRemoveVTopofBMS(50, 0);
                if (best_removed_v == -1){
                    if (undomPointArray->size() == 0){
                        UpdateBestSolution();
                        isShutDownSearchStep = true;
                        return;
                    }
                    else{
                        chooseMarkCuttree();
                        break;
                    }
                }
                else{
                    Remove(best_removed_v);
                    chooseMarkCuttree();
                    break;
                }
            }
            best_removed_v = ChooseRemoveVTopofBMS(50, 1);
            if (best_removed_v != -1){
                Remove(best_removed_v);
                removeUpdate(best_removed_v);
                time_stamp[best_removed_v] = step;
            }
        }
        bool flag_break = false;
        while (undomPointArray->size() != 0 && currentWeight < bestWeight){
            int best_add_v = ChooseAddVsubscorefast(50);
            if (best_add_v == -1){
                flag_break = true;
                break;
            }
            Add(best_add_v);
            addUpdate(best_add_v);
            time_stamp[best_add_v] = step;
        }
        if (currentWeight < bestWeight && flag_break == false){updateSolution();}
        if (undomPointArray->size() == 0 && currentWeight < bestWeight){
            neighborSize = 1;
            improvementCount++;
            NOimprovementstep = 0;
            TreeNoimprovementstep = 0;
        }
        else{
            for (size_t i = 0; i < undomPointArray->size(); i++){addWeight(undomPointArray->element_at(i));}
            totalweight += undomPointArray->size();
            NOimprovementstep++;
            if (neighborSize > maxNeighborSize){neighborSize = 1;}
            else{neighborSize++;}
        }
        TreeNoimprovementstep++;
        stepAction++;
        step++;
        if (NOimprovementstep % NoImpro_trystep == 0 && NOimprovementstep > 0){chooseMarkCuttree();}
    }
}


void Framework2TarjanFocus(){
    int neighborSize = tarjanRemoveSize;
    int currentLoop = 1;
    MarkCut();
    while (currentLoop < tarjanLoopCount * tarjanLoopCoefficient){
        if (candidateArray->size() == 1){return;}
        removedNodeNeighbor->clear();
        for (size_t i = 0; i < neighborSize; i++){
            int best_removed_v = -1;
            if (i == 0){best_removed_v = ChooseRemoveVTopofBMS(100, 0);}
            else{
                if (!removedNodeNeighbor->empty() && !candidateArray->empty()){best_removed_v = ChooseRemoveVFromArray(removedNodeNeighbor, 0);}
            }
            if (best_removed_v != -1){
                Remove(best_removed_v);
                if (!candidateArray->empty() && i < neighborSize - 1){MarkCut();}
                for (int n = 0; n < v_degree[best_removed_v]; ++n){
                    int neighbor = v_adj[best_removed_v][n];
                    if (v_in_c[neighbor] == 1 && !removedNodeNeighbor->is_in_array(neighbor)){removedNodeNeighbor->insert_element(neighbor);}
                    for (int m = 0; m < v_degree[neighbor]; m++){
                        int neinei = v_adj[neighbor][m];
                        if (v_in_c[neinei] == 1 && !removedNodeNeighbor->is_in_array(neinei)){removedNodeNeighbor->insert_element(neinei);}
                    }
                }
                time_stamp[best_removed_v] = step;
            }
        }
        while (undomPointArray->size() != 0 && currentWeight < bestWeight){
            int best_add_v = ChooseAddVsubscorefast(50, 0);
            if (best_add_v == -1 || time_stamp[best_add_v] == step)
                break;
            Add(best_add_v);
            RemoveRedundant();
            time_stamp[best_add_v] = step;
        }
        if (currentWeight < bestWeight){updateSolution();}
        if (undomPointArray->size() == 0 && currentWeight < bestWeight){
            NOimprovementstep = 0;
            Noimprovementstepfocus = 0;
        }
        else{
            for (size_t i = 0; i < undomPointArray->size(); i++){addWeight(undomPointArray->element_at(i));}
            totalweight += undomPointArray->size();
            NOimprovementstep++;
            Noimprovementstepfocus++;
        }
        step++;
        MarkCut();
        currentLoop++;
    }
}
#endif
