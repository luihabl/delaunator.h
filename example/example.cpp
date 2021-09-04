#include <delaunator.h>

#include <algorithm>
#include <set>
#include <utility>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>
#include <memory>
#include <string>
#include <sstream>

using namespace std; 

std::vector<double> read_csv(std::string filename)
{
    ifstream data(filename);
    string str;

    std::vector<double> d;
    while (getline(data, str))
    {
        double x, y;
        istringstream iss(str);
        string token;

        getline(iss, token, ',');
        x = std::stod(token);

        getline(iss, token, ',');
        y = std::stod(token);

        d.push_back(x);
        d.push_back(y);
    }

   return d;
}


class Graph{
private:
        int V,E;
        vector<pair<double, pair<int,int> > > edges;
        vector<pair<double, pair<int,int> > > MST;
public:
        Graph(int V,int E);
        void addEdge(int u,int v, double w);
        double kruskalMST();
        void printMST();
        vector<pair<double ,pair<int,int>>> getMST();
};
Graph::Graph(int V,int E){
    this->V = V;
    this->E = E;
}
void Graph::addEdge(int u,int v, double w){
    edges.push_back({w,{u,v}});
}

void Graph::printMST(){
    vector<pair<double ,pair<int,int> > >::iterator it;
    for(it = MST.begin();it!=MST.end();it++){
        cout << it->second.first << " - " << it->second.second << endl;
    }
}

vector<pair<double ,pair<int,int>>> Graph::getMST()
{
    return MST;
}


struct DisjointSet{
    int *parent,*rnk;
    int n;

    DisjointSet(int n){
        this->n = n;
        parent = new int[n+1];
        rnk = new int[n+1];

        for(int i=0;i<=n;i++){
            rnk[i] = 0;
            parent[i] = i;
        }
    }
    int Find(int u){
        if(u == parent[u]) return parent[u];
        else return Find(parent[u]);
    }

    void Union(int x,int y){
        x = Find(x);
        y = Find(y);
        if(x != y){
            if(rnk[x] < rnk[y]){
                rnk[y] += rnk[x];
                parent[x] = y;
            }
            else{
                rnk[x] += rnk[y];
                parent[y] = x;
            }
        }
    }
    ~DisjointSet() {
        delete parent;
        delete rnk;
    }
};

double Graph::kruskalMST(){
    double MSTWeight = 0; //sum of all vertex weights
    sort(edges.begin(),edges.end());
    //for all u in G_v
    //    MAKE-SET(u)
    DisjointSet ds(this->V);

    vector<pair<double,pair<int,int> > >::iterator it;
    // for all edges in G
    for(it = edges.begin(); it!=edges.end();it++){
        int u = it->second.first;
        int v = it->second.second;

        int setU = ds.Find(u);
        int setV = ds.Find(v);


        if(setU != setV){
            double w = it->first;
            MST.push_back({w,{u,v}});
            MSTWeight += it->first;

            ds.Union(setU,setV);
        }
    }
    return MSTWeight;
}

int main()
{   

    // std::vector<float> points = {
    //     -12, 8,
    //     -10, 6,
    //     -12, 4,
    //     12, 8,
    //     10, 6,
    //     12, 4
    // };
    // std::cout << "A" << std::endl;

    auto coords = read_csv("<path>");
    size_t np = coords.size() / 2;
    
    // std::cout << "B" << std::endl;
    del2d::Delaunator d(coords);

    // std::cout << d.triangles.size() << std::endl;

    // std::cout << "C" << std::endl;

    std::set<std::pair<size_t, size_t>> s;

    for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
        size_t a = d.triangles[i];
        size_t b = d.triangles[i+1];
        size_t c = d.triangles[i+2];

        s.insert(std::pair<size_t, size_t>(std::min(a,b), std::max(a,b)));
        s.insert(std::pair<size_t, size_t>(std::min(b,c), std::max(b,c)));
        s.insert(std::pair<size_t, size_t>(std::min(a,c), std::max(a,c)));       
    }

    std::vector<double> distances;

    const auto& dist = [](double x0, double y0, double x1, double y1)
    {
        return sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
    };
    

    for (auto const& el : s)
    {
        distances.push_back(dist(coords[el.first * 2 + 0], coords[el.first * 2 + 1], coords[el.second * 2 + 0], coords[el.second * 2 + 1]));
    }
    
    size_t i = 0;
    // for (auto const& el : s)
    // {
    //     std::cout << el.first << " " << el.second << " : " << distances[i] << std::endl;
    //     i++;
    // } 

    // std::cout << "D" << std::endl;


    Graph g(np, distances.size());

    i = 0;
    for (auto const& el : s)
    {
        g.addEdge(el.first, el.second, distances[i]);
        i++;
    } 

    double weight = g.kruskalMST();
    // g.printMST();

    auto mst = g.getMST();
    
    double m = 0;
    for(const auto& n : mst)
    {   
        m = std::max(m, n.first);
    }

    std::cout << "Result: " << m << std::endl;

}