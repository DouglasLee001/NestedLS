#include "basic.h"
#include <queue>
/*******update*******/

void ResetCandidate_small(){
    int j = 0;
    for (int v = 1; v < v_num + 1; v++){
        if (v_in_c[v] == 1 && v_fixed[v] == 0) {
            candidate[j] = v;
            index_in_candidate[v] = j;
            j++;
        }
        else{
            index_in_candidate[v] = 0;
            if (dominated[v] > 0){greyPointArray->insert_element(v);}
        }
    }
    candidate_size = j;
}

void UpdateBestSolution_small(){
    if ((int)currentWeight < (int)bestWeight){
        best_c_size = c_size;
        best_comp_time = TimeElapsed();
        best_step = step;
        bestWeight = currentWeight;
        memcpy(best_v_in_c, v_in_c, sizeof(int) * (v_num + 1));
//         cout << "weight: " << currentWeight  << "\ttime: " << best_comp_time << endl;
    }
}

void cutPointNoRecur_small(int root){
    int k, u, v, inde, top = 0;
    inde = 0;
    for (u = 0; u < candidate_size; u++){SF[candidate[u]] = first[candidate[u]];}
    for (u = 0; u < fixedNum; u++){SF[fixedSet[u]] = first[fixedSet[u]];}
    Stack[top] = root;
    dnf[root] = low[root] = ++inde;
    while (top >= 0){
        k = SF[Stack[top]];
        if (k != -1){
            if (v_in_c[edge[k].v] == 0){
                SF[Stack[top]] = edge[k].next;
                continue;
            }
            SF[Stack[top]] = edge[k].next;
            if (dnf[edge[k].v] == 0){
                child[Stack[top]]++;
                Stack[++top] = edge[k].v;
                low[edge[k].v] = dnf[edge[k].v] = ++inde;
            }
            else{low[Stack[top]] = (low[Stack[top]] < dnf[edge[k].v]) ? low[Stack[top]] : dnf[edge[k].v];}
        }
        else{
            if (top > 0){
                u = Stack[top - 1];
                v = Stack[top];
                low[u] = (low[u] < low[v]) ? low[u] : low[v];
                if ((u != root && low[v] >= dnf[u]) || (u == root && child[u] >= 2)){
                    isCut[u] = 1;
                    cutPointSet[cutIndex++] = u;
                }
            }
            top--;
        }
    }
}
void minusWeight_small(int node, int tobeminus);
void updateWeight_small(){
    double avg_weight = totalweight / v_num;
    for (int i = 1; i <= v_num; i++){
        int newweight = weight_para_aphle * frequency[i] + avg_weight * (1 - weight_para_aphle);
        if (newweight < 1)
            newweight = 1;
        int tobeminus = (frequency[i] - newweight);
        minusWeight_small(i, tobeminus);
        totalweight -= tobeminus;
    }
}

/******construct******/
void resetScore_small(){
    fill_n(score, v_num + 1, 0);
    fill_n(subscore, v_num + 1, 0);
    for (size_t i = 1; i < v_num + 1; i++){
        score[i] = v_degree[i] + 1;
        subscore[i] = frequency[i];
        for (size_t j = 0; j < v_degree[i]; j++){subscore[i] += frequency[v_adj[i][j]];}
    }
}

void constructIncreaseDominate_small(int v, int v_dominater){
    if (dominated[v] == 0){
        --score[v];
        subscore[v] -= frequency[v];
        for (size_t i = 0; i < v_degree[v]; i++){
            int u = v_adj[v][i];
            --score[u];
            subscore[u] -= frequency[v];
        }
        Dom(v);
        onlydominate[v] = v_dominater;
        if (v != v_dominater){greyPointArray->insert_element(v);}
    }
    else if (dominated[v] == 1){
        if (v_in_c[v] == 1){
            ++score[v];
            subscore[v] += frequency[v];
        }
        else{
            int v_dominater = onlydominate[v];
            ++score[v_dominater];
            subscore[v_dominater] += frequency[v];
        }
    }
    ++dominated[v];
}

void constructAdd_small(int v, bool isInitial = false){
    if (v_in_c[v] == 1){return;}
    v_in_c[v] = 1;
    if (isInitial){last_v_in_c[v] = 1;}
    c_size++;
    currentWeight += weight[v];
    greyPointArray->delete_element(v);
    int new_score = -score[v];
    int new_subscore = -subscore[v];
    constructIncreaseDominate_small(v, v);
    for (size_t i = 0; i < v_degree[v]; i++){
        int u = v_adj[v][i];
        constructIncreaseDominate_small(u, v);
    }
    score[v] = new_score;
    subscore[v] = new_subscore;
}

int NewSolutionChooseVFromMethodA_small()
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
            else if (cscore == best_cscore){
                if (subWeight[i] > subWeight[best_add_v]){best_add_v = i;}
            }
        }
    }
    else{
        for (size_t i = 0; i < greyPointArray->size(); i++){
            int greyPoint = greyPointArray->element_at(i);
            cscore = (int)(ave_weight * subscore[greyPoint] / weight[greyPoint]);
            if (cscore > best_cscore){
                best_add_v = greyPoint;
                best_cscore = cscore;
            }
            else if (cscore == best_cscore){
                if (subWeight[greyPoint] > subWeight[best_add_v]){best_add_v = greyPoint;}
            }
        }
    }
    return best_add_v;
}

void ConstructRestartSolution_small(){
    undomPointArray->clear();
    for (int v = 1; v < v_num + 1; ++v){undomPointArray->insert_element(v);}
    resetScore_small();
    while (undomPointArray->size() != 0){
        int addV = -1;
        addV = NewSolutionChooseVFromMethodA_small();
        if (addV != -1){
            constructAdd_small(addV);
        }
    }
}

double getNewSolutionRepeatDegree_small(){
    int lastCount = 0;
    int newCount = 0;
    int repeatCount = 0;
    for (size_t v = 1; v < v_num + 1; v++){
        if (v_in_c[v] == 1 && last_v_in_c[v] == 1){repeatCount++;}
        else if (v_in_c[v] == 1){newCount++;}
        else if (last_v_in_c[v] == 1){lastCount++;}
    }
    return repeatCount * 1.0 / (lastCount + repeatCount + newCount) * 100;
}

void updateRedundantV_small(int v){
    if (!moduleRemoveRedundant){return;}
    if (v_in_c[v] == 1 && score[v] == 0 && v_fixed[v] == 0){redundantNodes->insert_element(v);}
    else{redundantNodes->delete_element(v);}
}
bool Remove_small(int v);
void MarkCut_small();
void RemoveRedundant_small(){
    if (!moduleRemoveRedundant){return;}
    MarkCut_small();
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
            Remove_small(bestRemoveV);
            MarkCut_small();
        }
        else break;
    }
}

void prepareForLocalSearch_small(){
    ConstructRestartSolution_small();
    bestWeightInTurn = currentWeight;
    candidate_size = 0;
    for (int i = 1; i <= v_num; ++i){candidate[candidate_size++] = i;}
    MarkCut_small();
    for (int i = 0; i < cutIndex; ++i){
        if (isCut[cutPointSet[i]] != 0){
            int cutPoint = cutPointSet[i];
            v_fixed[cutPoint] = 1;
            fixedSet[fixedNum++] = cutPoint;
        }
    }
    ResetCandidate_small();
    RemoveRedundant_small();
    UpdateBestSolution_small();
}

void Restart_small(){
    stepAction = 1;
    NOimprovementstep = 0;
    Noimprovementstepfocus = 0;
    greyPointArray->clear();
    fill_n(v_in_c, v_num + 1, 0);
    fill_n(dominated, v_num + 1, 0);
    fill_n(taburemove, v_num + 1, 0);
    fill_n(tabuadd, v_num + 1, 0);
    fill_n(conf_change, v_num + 1, 1);
    fill_n(v_fixed, v_num + 1, 0);
    fill_n(isCut, v_num + 1, 0);
    fill_n(cutPointSet, v_num + 1, 0);
    c_size = 0;
    cutIndex = 0;
    currentWeight = 0;
    if (lastRepeatDegree > restartUpdateFrequencyThreshold){updateWeight_small();}
    ConstructRestartSolution_small();
    fill_n(subWeight, v_num + 1, 0);
    ResetCandidate_small();
    RemoveRedundant_small();
    lastRepeatDegree = getNewSolutionRepeatDegree_small();
    memcpy(last_v_in_c, v_in_c, sizeof(int) * (v_num + 1));
    bestWeightInTurn = currentWeight;
    MarkCut_small();
}

void increase_dominate_small(long v, long source_v){
    if (dominated[v] == 0){
        --score[v];
        subscore[v] -= frequency[v];
        updateRedundantV_small(v);
        for (int i = 0; i < v_degree[v]; ++i){
            long u = v_adj[v][i];
            --score[u];
            subscore[u] -= frequency[v];
            updateRedundantV_small(v);
        }
        Dom(v);
        onlydominate[v] = source_v;
        greyPointArray->insert_element(v);
    }
    else if (dominated[v] == 1){
        if (v_in_c[v] == 1){
            ++score[v];
            subscore[v] += frequency[v];
            updateRedundantV_small(v);
        }
        else{
            long u = onlydominate[v];
            ++score[u];
            subscore[u] += frequency[v];
            updateRedundantV_small(u);
        }
    }
    ++dominated[v];
}

bool Add_small(int v){
    if (v_in_c[v] == 1 || v < 0)
        return false;
    greyPointArray->delete_element(v);
    if (currentMode == ChooseMode::ModeA){subWeight[v]++;}
    int new_score = -score[v];
    int new_subscore = -subscore[v];
    increase_dominate_small(v, v);
    for (int i = 0; i < v_degree[v]; ++i){
        long u = v_adj[v][i];
        increase_dominate_small(u, v);
    }
    score[v] = new_score;
    subscore[v] = new_subscore;
    v_in_c[v] = 1;
    updateRedundantV_small(v);
    currentWeight += weight[v];
    c_size++;
    candidate[candidate_size] = v;
    index_in_candidate[v] = candidate_size++;
    int u;
    int tmp = conf_change[v];
    for (int i = 0; i < v_degree[v]; i++){
        u = v_adj[v][i];
        for (int j = 0; j < v_degree[u]; j++)
            conf_change[v_adj[u][j]] = 1;
        conf_change[u] = 1;
    }
    conf_change[v] = tmp;
    taburemove[v] = (step + 1 + rand() % 3);
    return true;

}

void decrease_dominate_small(int v){
    if (dominated[v] == 1){
        ++score[v];
        subscore[v] += frequency[v];
        updateRedundantV_small(v);
        for (int i = 0; i < v_degree[v]; ++i){
            int u = v_adj[v][i];
            ++score[u];
            subscore[u] += frequency[v];
            updateRedundantV_small(u);
        }
        Undom(v);
        greyPointArray->delete_element(v);
    }
    else if (dominated[v] == 2){
        if (v_in_c[v]){
            --score[v];
            subscore[v] -= frequency[v];
            updateRedundantV_small(v);
        }
        else{
            for (int i = 0; i < v_degree[v]; ++i){
                long u = v_adj[v][i];
                if (v_in_c[u]){
                    --score[u];
                    subscore[u] -= frequency[v];
                    updateRedundantV_small(u);
                    onlydominate[v] = u;
                    break;
                }
            }
        }
    }
    --dominated[v];
}

bool Remove_small(int v){
    greyPointArray->insert_element(v);
    v_in_c[v] = 0;
    currentWeight -= weight[v];
    c_size--;
    int last_candidate_v = candidate[--candidate_size];
    int index = index_in_candidate[v];
    candidate[index] = last_candidate_v;
    index_in_candidate[last_candidate_v] = index;
    int new_score = -score[v];
    int new_subscore = -subscore[v];
    decrease_dominate_small(v);
    int neighbours = v_degree[v];
    for (int i = 0; i < neighbours; ++i){
        int u = v_adj[v][i];
        decrease_dominate_small(u);
    }
    score[v] = new_score;
    subscore[v] = new_subscore;
    updateRedundantV_small(v);
    int u;
    for (int i = 0; i < v_degree[v]; i++){
        u = v_adj[v][i];
        for (int j = 0; j < v_degree[u]; j++)
            conf_change[v_adj[u][j]] = 1;
        conf_change[u] = 1;
    }
    conf_change[v] = 0;
    tabuadd[v] = step + 1;
    return true;
}


void addWeight_small(int node){
    int increment = 5;
    frequency[node] += increment;
    for (int i = 0; i < v_degree[node]; i++){subscore[v_adj[node][i]] += increment;}
    subscore[node] += increment;
}

void minusWeight_small(int node, int tobeminus){
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

bool cmp_small(int a, int b){return (subscore[a] / weight[a] > subscore[b] / weight[b]);}

int ChooseRemoveVTopofBMS_small(int count, int choice){
    int v, i;
    int best_score = -weightthreshold;
    int best_remove_v = -1;
    const int toberemovedNum1 = (int)(count);
    int toberemoved1[toberemovedNum1];
    int topIndex = 0;
    int cscore;
    for (i = 0; i < count; i++){
        v = candidate[rand() % candidate_size];
        if (v_fixed[v] == 1 || isCut[v] == 1){continue;}
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
        if (topIndex > 0){
            sort(toberemoved1, toberemoved1 + topIndex, cmp_small);
            int topof = (topIndex * 0.4) + 1;
            return toberemoved1[rand() % topof];
        }
        else{return -1;}
    }
}

int ChooseRemoveVFromArray_small(Array *removedNodeNeighbor, int choice){ //choice==0,tarjan;choice==1,tree
    int best_remove_v = -1;
    int cscore;
    int best_cscore = -weightthreshold;
    for (size_t i = 0; i < removedNodeNeighbor->size(); i++){
        int v = removedNodeNeighbor->element_at(i);
    if (v_in_c[v] == 0 || v_fixed[v] == 1 || isCut[v] == 1){continue;}
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

int ChooseAddVsubscorefast_small(int choice=0){
    int base_v, add_v;
    int cscore;
    int best_score = -weightthreshold;
    int best_add_v = -1;
    map<int, int> m;
    const int tobeaddNum1 = (int)(undomPointArray->size() * v_degree[maxDegreeNode]);
    int tobeadd1[tobeaddNum1];
    int topIndex = 0;
    for (int i = 0; i < undomPointArray->size(); ++i){
        base_v = undomPointArray->element_at(i);
        for (int j = 0; j < v_degree[base_v]; ++j){
            add_v = v_adj[base_v][j];
            if (m.find(add_v) == m.end() && greyPointArray->is_in_array(add_v) && (weight[add_v] + currentWeight < bestWeightInTurn)){
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
            sort(tobeadd1, tobeadd1 + topIndex, cmp_small);
            int topof = (topIndex * 0.4) + 1;
            return tobeadd1[rand() % topof];
        }
        else{return -1;}
    }
}

void MarkCut_small(){
    ind = 0;
    for (int i = 0; i < cutIndex; i++){isCut[cutPointSet[i]] = 0;}
    cutIndex = 0;
    if (candidate_size == 0){return;}
    root = candidate[0];
    cutPointNoRecur_small(root);
    for (int i = 0; i < candidate_size; i++){child[candidate[i]] = low[candidate[i]] = dnf[candidate[i]] = 0;}
    for (int i = 0; i < fixedNum; i++){child[fixedSet[i]] = low[fixedSet[i]] = dnf[fixedSet[i]] = 0;}
}

void Framework2TarjanScatter_small(){
    try_step = 1000;
    int improvementCount = 0;
    int neighborSize = 1;
    MarkCut_small();
    while (true){
        if (stepAction % try_step == 0){
            int timenow = TimeElapsed();
            if (timenow > cutoff_time){return;}
        }
        if (stepAction >= (100 * floor0)){return;}
        if (NOimprovementstep >= (floor0 * 40)){return;}
        if (c_size > 50 && Noimprovementstepfocus >= (floor0 * 10) * betweemNoimpr){return;}
        if (candidate_size == 1){return;}
        for (size_t i = 0; i < neighborSize; i++){
            int best_removed_v = -1;
            best_removed_v = ChooseRemoveVTopofBMS_small(100, 0);
            if (best_removed_v != -1){
                Remove_small(best_removed_v);
                if (candidate_size != 0 && i < neighborSize - 1){MarkCut_small();}
                time_stamp[best_removed_v] = step;
            }
        }
        while (undomPointArray->size() != 0 && currentWeight < bestWeightInTurn){
            int best_add_v = ChooseAddVsubscorefast_small(0);
            if (best_add_v == -1 || time_stamp[best_add_v] == step){break;}
            Add_small(best_add_v);
            RemoveRedundant_small();
            time_stamp[best_add_v] = step;
        }
        if (undomPointArray->size() == 0 && currentWeight < bestWeightInTurn){
            bestWeightInTurn = currentWeight;
            UpdateBestSolution_small();
            neighborSize = 1;
            improvementCount++;
            NOimprovementstep = 0;
            Noimprovementstepfocus = 0;
            if (c_size > 50 && betweemNoimpr > 1){betweemNoimpr /= 2;}
        }
        else{
            for (size_t i = 0; i < undomPointArray->size(); i++){addWeight_small(undomPointArray->element_at(i));}
            totalweight += undomPointArray->size();
            NOimprovementstep++;
            Noimprovementstepfocus++;
            if (neighborSize > maxNeighborSize || neighborSize >= candidate_size - 1){neighborSize = 1;}
            else{neighborSize++;}
        }
        stepAction++;
        step++;
        MarkCut_small();
    }
}

void localSearchFramework2_small(){
    while (TimeElapsed() < cutoff_time){
        Framework2TarjanScatter_small();
        if (NOimprovementstep >= (40 * floor0) || (stepAction >= (100 * floor0)) || candidate_size == 1){Restart_small();}
        if (c_size > 50 && Noimprovementstepfocus >= betweemNoimpr * (10 * floor0)){
            Noimprovementstepfocus = 0;
            if (betweemNoimpr < 4){betweemNoimpr *= 2;}
        }
    }
}

void enter_ls(){
    if(v_num<=1000){
        moduleRemoveRedundant=true;
        prepareForLocalSearch_small();
        localSearchFramework2_small();
    }
    else{
        prepareForLocalSearch();
        localSearchFramework2();
    }
}
