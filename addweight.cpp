#include <iostream>
#include <cstdio>
#include <stack>
#include <cstring>
#include <queue>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <random>
using namespace std;
using std::default_random_engine;
using std::uniform_int_distribution;

struct edgeNode
{
    int v1, v2;
};
edgeNode *edges;
long **adj;
long *visit;
long *degree;
vector<long> nowCnt;
vector<long> bigCnt;
// vector<int> weight;
long v, e;
long edge = 0;
long ind = 0;

void freeMemory()
{
    for (int i = 0; i < v + 1; ++i)
    {
        delete[] adj[i];
    }
    delete[] adj;
    delete[] visit;
    delete[] degree;
    delete[] edges;
}

void bfs(long root)
{
    nowCnt = {};
    queue<long> q;
    nowCnt.push_back(root);
    q.push(root);
    visit[root] = 1;
    ind++;
    long cur;
    while (!q.empty())
    {
        cur = q.front();
        q.pop();
        for (int i = 0; i < degree[cur]; i++)
        {
            long tt = adj[cur][i];
            if (visit[tt] == 0)
            {
                visit[tt] = 1;
                q.push(tt);
                nowCnt.push_back(tt);
                ind++;
            }
        }
    }
}

int main(int argc, char *argv[])
{
    string filename = argv[1];
    freopen(argv[1], "r", stdin);
    long a, b;
    string tmp1, tmp2;
    cin >> tmp1 >> tmp2 >> v >> e;

    visit = new long[v + 1];
    degree = new long[v + 1];
    adj = new long *[v + 1];
    edges = new edgeNode[e];
    memset(degree, 0, sizeof(long) * (v + 1));
    memset(visit, 0, sizeof(long) * (v + 1));
    

    for (long i = 0; i < e; i++)
    {
        cin >> tmp1 >> a >> b;
        edges[i].v1 = a;
        edges[i].v2 = b;
        degree[a]++;
        degree[b]++;
    }
    adj[0] = nullptr;
    for (int i = 1; i < v + 1; ++i)
    {
        adj[i] = new long[degree[i]];
    }
    int *degree_tmp = new int[v + 1];
    memset(degree_tmp, 0, v + 1);
    for (int i = 0; i < e; ++i)
    {
        int v1 = edges[i].v1;
        int v2 = edges[i].v2;
        adj[v1][degree_tmp[v1]++] = v2;
        adj[v2][degree_tmp[v2]++] = v1;
    }
    delete[] degree_tmp;
    fclose(stdin);
    // find max cnt
    for (int i = 1; i < v + 1; ++i)
    {
        if (degree[i] > 0 && visit[i] == 0)
        {
            bfs(i);
            if (nowCnt.size() > bigCnt.size())
            {
                bigCnt.assign(nowCnt.begin(), nowCnt.end());
            }
        }
    }
    edge = 0;
    map<long, long> m;
    sort(bigCnt.begin(), bigCnt.end());
    for (size_t i = 0; i < bigCnt.size(); i++)
    {
        long key = bigCnt[i];
        long value = i + 1;
        m.insert(std::pair<long, long>(key, value));
        for (size_t j = 0; j < degree[bigCnt[i]]; j++)
        {
            if (adj[bigCnt[i]][j] > bigCnt[i])
            {
                edge++;
            }
        }
    }

    freopen(argv[2], "w", stdout);
    cout << "p edge " << bigCnt.size() << ' ' << edge << endl;
    for (long i : bigCnt)
    {
        cout << "v " << m[i] << " " << m[i]%200+1 << endl;
    }

    for (long &i : bigCnt)
    {
        for (int j = 0; j < degree[i]; j++)
            if (adj[i][j] > i)
                cout << "e " << m[i] << ' ' << m[adj[i][j]] << endl;
    }

    fclose(stdout);
    freeMemory();
    return 0;
}
