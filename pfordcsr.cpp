#include <mpi.h>
#include <iostream>
#include <vector>
#include <queue>
#include <cstring>
#include <fstream>

using namespace std;

#define V 1632405

struct CSRGraph {
    vector<int> vertices;  
    vector<int> edges;     
};

bool bfs(CSRGraph& graph, int s, int t, vector<int>& parent, int rank) {
    bool* visited = new bool[V];
    memset(visited, 0, sizeof(bool) * V);

    queue<int> q;
    int path_flow = 0;
    if (rank == 0) {
        q.push(s);
        visited[s] = true;
        parent[s] = -1;
    }

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (int i = graph.vertices[u]; i < graph.vertices[u + 1]; i++) {
            int v = graph.edges[i];
            if (visited[v] == false) {
                if (v == t) {
                    parent[v] = u;
                    delete[] visited;
                    return true;
                }
                if (rank == 0) {
                    q.push(v);
                    parent[v] = u;
                    visited[v] = true;
                    path_flow = 1;
                    cout << "MPI Send " << rank << endl;
                    MPI_Send(visited, V, MPI_C_BOOL, 1, 0, MPI_COMM_WORLD);
                    cout << "MPI Send " << rank << endl;
                    MPI_Send(&path_flow, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
                } else {
                    cout << "MPI " << rank << endl;
                    MPI_Recv(visited, V, MPI_C_BOOL, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    cout << "MPI" << rank << endl;
                    MPI_Recv(&path_flow, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }
    }

    delete[] visited;
    return false;
}

int fordFulkerson(CSRGraph& graph, int s, int t, int rank) {
    vector<int> parent(V, -1);
    int max_flow = 0;

    while (bfs(graph, s, t, parent, rank)) {
        int path_flow = 1;

        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            if (rank == 0) {

                path_flow = 1;
                MPI_Send(&path_flow, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
            } else {

                MPI_Recv(&path_flow, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        max_flow += path_flow;
    }

    return max_flow;
}

int main(int argc, char* argv[]) {
    printf("#Hi hello\n");
    printf("#Debug 2\n");

    MPI_Init(NULL, NULL);
    printf("#DEBUG 1.5\n");

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("DEBUG 2 After RANK AND SIZE\n");

    if (argc != 2) {
        if (rank == 0) {
            cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    string inputFile = argv[1];
    ifstream file(inputFile);

    int source = 0;
    int sink = 1632404;

    CSRGraph graph;

    graph.vertices.push_back(0); 

    int from, to;

    while (file >> from >> to) {

        graph.edges.push_back(to);
    }


    for (int i = 1; i <= V; i++) {
        int row_start = graph.vertices[i - 1] + (i - 1);
        graph.vertices.push_back(row_start);
    }

    if (rank == 0) {
        int max_flow = fordFulkerson(graph, source, sink, rank);
        cout << "The maximum possible flow is " << max_flow << endl;
    } else {
        fordFulkerson(graph, source, sink, rank);
    }

    MPI_Finalize();

    return 0;
}
