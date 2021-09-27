#ifndef _CHECKER_H
#define _CHECKER_H

#include "basic.h"

bool test_score()
{
    vector<int> correct_score(v_num + 1, 0);
    vector<int> correct_dominate(v_num + 1, 0);
    for (int v = 1; v < v_num + 1; ++v)
    {
        if (v_in_c[v])
        {
            correct_dominate[v]++;
            for (int i = 0; i < v_degree[v]; ++i)
            {
                int u = v_adj[v][i];
                correct_dominate[u]++;
            }
        }
    }
    for (int v = 1; v < v_num + 1; ++v)
    {
        if (correct_dominate[v] == 0)
        {
            correct_score[v]++;
            for (int i = 0; i < v_degree[v]; ++i)
            {
                int u = v_adj[v][i];
                correct_score[u]++;
            }
        }
        else if (correct_dominate[v] == 1)
        {
            if (v_in_c[v])
            {
                correct_score[v]--;
            }
            else
            {
                for (int i = 0; i < v_degree[v]; ++i)
                {
                    int u = v_adj[v][i];
                    if (v_in_c[u])
                    {
                        correct_score[u]--;
                    }
                }
            }
        }
    }
    bool right = true;
    for (int v = 1; v < v_num + 1; ++v)
    {
        if (correct_dominate[v] != dominated[v])
        {
            std::cout << "dominated wrong: " << v << std::endl;
            std::cout << "correct: " << correct_dominate[v]
                      << ", wrong: " << dominated[v] << std::endl;
            right = false;
        }
        if (correct_score[v] != score[v])
        {
            std::cout << "score wrong: " << v << std::endl;
            std::cout << "correct: " << correct_score[v] << ", wrong: " << score[v]
                      << std::endl;
            std::cout << v_in_c[v] << std::endl;
            right = false;
        }
    }
    if (!right)
        abort();

    right = true;
    int undom_num = 0;
    for (int v = 1; v < v_num + 1; ++v)
    {
        if (dominated[v] == 0)
        {
            undom_num++;
            bool okay = false;
            for (int i = 0; i < undomPointArray->size(); ++i)
            {
                if (undomPointArray->element_at(i) == v)
                {
                    okay = true;
                    break;
                }
            }
            if (!okay)
            {
                std::cout << "not in undom_stack: " << v << std::endl;
                right = false;
            }
        }
    }
    if (undom_num != undomPointArray->size())
    {
        std::cout << "undom_num wrong!" << std::endl;
    }
    if (!right)
        abort();

    return true;
}
bool CheckSolution()
{
    int v;
    memset(dominated, 0, sizeof(int) * (v_num + 1));

    for (v = 1; v < v_num + 1; v++)
    {
        if (best_v_in_c[v])
            for (int j = 0; j < v_degree[v]; j++)
                dominated[v_adj[v][j]]++;
        dominated[v]++;
    }
    for (v = 1; v < v_num + 1; v++)
        if (dominated[v] == 0)
            return false;
    return true;
}

// check if the solution is connected
bool CheckSolutionIsConnected()
{
    queue<int> q;
    vector<int> visit(v_num + 1, 0);
    for (int v = 1; v < v_num + 1; ++v)
    {
        if (best_v_in_c[v])
        {
            q.push(v);
            visit[v]++;
            break;
        }
    }

    int visitCount = 0;
    while (!q.empty())
    {
        int v = q.front();
        q.pop();
        visitCount++;
        for (int n = 0; n < v_degree[v]; ++n)
        {
            int neighbor = v_adj[v][n];
            if (best_v_in_c[neighbor] && visit[neighbor] == 0)
            {
                q.push(neighbor);
                visit[neighbor]++;
            }
        }
    }
    return visitCount == best_c_size;
}

bool CheckGraphIsConnected()
{
    queue<int> q;
    vector<int> visit(v_num + 1, 0);
    q.push(1);
    visit[1] = 1;
    int visitCount = 0;
    while (!q.empty())
    {
        int v = q.front();
        q.pop();
        visitCount++;
        for (int n = 0; n < v_degree[v]; ++n)
        {
            int neighbor = v_adj[v][n];
            if (visit[neighbor] == 0)
            {
                q.push(neighbor);
                visit[neighbor]++;
            }
        }
    }
    return visitCount == v_num;
}

bool checkCurrentSolutionIsconnected()
{
    queue<int> q;
    vector<int> visit(v_num + 1, 0);
    for (int v = 1; v < v_num + 1; ++v)
    {
        if (v_in_c[v])
        {
            q.push(v);
            visit[v]++;
            break;
        }
    }

    int visitCount = 0;
    while (!q.empty())
    {
        int v = q.front();
        q.pop();
        visitCount++;
        for (int n = 0; n < v_degree[v]; ++n)
        {
            int neighbor = v_adj[v][n];
            if (v_in_c[neighbor] && visit[neighbor] == 0)
            {
                q.push(neighbor);
                visit[neighbor]++;
            }
        }
    }
    return visitCount == c_size;
}

#endif