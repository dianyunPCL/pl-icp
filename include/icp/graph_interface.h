#ifndef GRAPH_INTERFACE_H
#define GRAPH_INTERFACE_H
#include <bits/stdc++.h>
#include <Eigen/Core>

namespace GRAPH {

    constexpr auto NUMBER_OF_POINTS = 1081;

    class Vertex
    {
    public:
        Vertex() = delete;
        Vertex(double tx, double ty, double rotate);
        Vertex(double *p);
        Vertex(const Vertex &other);
        Vertex(Vertex &&other) = default;
        Vertex &operator=(const Vertex &other);
        Vertex &operator=(Vertex &&other) = default;
        ~Vertex() = default;

        double position[3];
    };

    class LaserData
    {
    public:
        LaserData() = delete;
        LaserData(int id, int *data, char *flag, double yaw);
        LaserData(const LaserData &other);
        LaserData(LaserData &&other) = default;
        ~LaserData();
        LaserData &operator=(const LaserData &other);
        LaserData &operator=(LaserData &&other) = default;
        void releaseData();

        int frame_id;
        int *data;
        char *flag;
        double yaw;
        double relativeTransfer[3];
        double globalPosition[3];
    };

    class Edge
    {
    public:
        Edge() = delete;
        Edge(int id_to, int id_from, double tx, double ty, double rotate);
        Edge(int id_to, int id_from, double *p);
        Edge(const Edge &other);
        Edge(Edge &&other) = default;
        Edge &operator=(const Edge &other);
        Edge &operator=(Edge &&other) = default;
        ~Edge() = default;

        int id_to;
        int id_from;
        double transfer[3];
    };

    template<typename V=Vertex, typename E=Edge, typename LD=LaserData>
    class GraphInterface
    {
    public:
        using vertexs = std::vector<V>;
        using unique_vertex = std::unique_ptr<vertexs>;
        using edges = std::vector<E>;
        using unique_edge = std::unique_ptr<edges>;
        using laserdatas = std::vector<LD>;
        using laserdata_iterator = typename laserdatas::iterator;

        virtual ~GraphInterface() = default;
        GraphInterface(const GraphInterface<V, E, LD> &other) : graph_vertexs(), graph_edges() {};
        virtual void add_vertex(int frame_id) = 0;
        virtual void add_edge(int frame_id_from, int frame_id_to, double rotate) = 0;
        virtual int add_laserdata(int *data, int *flag, double yaw);
        void updatePosition(laserdata_iterator begin);
        decltype(auto) beginToOptimize(unique_vertex &pVertex, unique_edge &pEdge);

        unique_vertex graph_vertexs;
        unique_edge graph_edges;
        laserdatas key_frames;
        int lastOpti;
        
    protected:
        GraphInterface() : graph_vertexs(new vertexs), graph_edges(new edges), lastOpti(0) {};
    };

}

#endif // !GRAPH_INTERFACE_H
