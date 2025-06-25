#pragma once
#include <Eigen/Core>

// Forward‐declares
class RootFeatures;
class Edge;
class Vertex;

// A generic callbacks namespace, templated by (FeatureType, EdgeType, VertexType)
namespace RAGCallbacks
{

    template <typename Feature, typename EdgeT, typename VertexT>
    using GetFeatures = Feature *(*)(int *, double *, int, int, int);

    template <typename Feature, typename EdgeT, typename VertexT>
    using UpdateFeatures = void (*)(Feature *, Feature *, EdgeT *);

    template <typename Feature, typename EdgeT, typename VertexT>
    using MergingCriteria = double (*)(Feature *, Feature *, EdgeT *);
}