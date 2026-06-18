#include "gi_discrete_gradient_computer.h"
#include "gi_morse_smale_complex_basic.h"

#include <cmath>
#include <iostream>
#include <map>
#include <random>
#include <utility>
#include <vector>

typedef GInt::MorseSmaleComplexBasic<float, GInt::Accurate2D::MeshType, GInt::Accurate2D::MeshFuncType, GInt::Accurate2D::GradType> MyMscType;

static std::map<std::pair<INDEX_TYPE, INDEX_TYPE>, int> endpointHistogram(MyMscType& msc) {
    std::map<std::pair<INDEX_TYPE, INDEX_TYPE>, int> hist;
    for (INT_TYPE aid = 0; aid < msc.numArcs(); ++aid) {
        const GInt::arc<float>& a = msc.getArc(aid);
        const INDEX_TYPE lowerCell = msc.getNode(a.lower).cellindex;
        const INDEX_TYPE upperCell = msc.getNode(a.upper).cellindex;
        ++hist[std::make_pair(lowerCell, upperCell)];
    }
    return hist;
}

int main() {
    const int rows = 1024;
    const int cols = 1024;

    std::vector<float> field(static_cast<size_t>(rows) * static_cast<size_t>(cols), 0.0f);
    std::mt19937 rng(123456789u);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    const float cx = 0.5f * static_cast<float>(cols);
    const float cy = 0.5f * static_cast<float>(rows);
    const float radius = 400.0f;
    const float r2 = radius * radius;
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            const float dx = static_cast<float>(c) - cx;
            const float dy = static_cast<float>(r) - cy;
            float v = dist(rng);
            if (dx * dx + dy * dy > r2) {
                v = 0.0f;
            }
            field[static_cast<size_t>(r) * static_cast<size_t>(cols) + static_cast<size_t>(c)] = v;
        }
    }

    GInt::Accurate2D::DiscreteGradientBuilder dgb;
    dgb.SetFloadArrayAndDims(cols, rows, field.data());
    dgb.SetNeededAccuracy(true, true);
    dgb.SetParallelism(4);
    dgb.ComputeDiscreteGradient();

    auto* mesh = dgb.GetTopoMesh();
    auto* meshfunc = dgb.GetMeshFunc();
    auto* grad = dgb.GetGrad();
    auto* gridfunc = dgb.GetGridFunc();

    const float maxval = gridfunc->GetMaxValue();
    const float minval = gridfunc->GetMinValue();
    const float persLimit = 0.1f * (maxval - minval);

    std::cout << "=== Generic no-geom ComputeFromGrad ===" << std::endl;
    MyMscType genericMSC(grad, mesh, meshfunc);
    genericMSC.SetBuildArcGeometry(GInt::Vec3b(false, false, false));
    genericMSC.ComputeFromGrad();
    const auto genericPreHist = endpointHistogram(genericMSC);
    const INT_TYPE genericNodes = genericMSC.numNodes();
    const INT_TYPE genericArcs = genericMSC.numArcs();
    genericMSC.ComputeHierarchy(persLimit);

    std::cout << "=== Dense-cache 2D no-geom builder ===" << std::endl;
    MyMscType denseMSC(grad, mesh, meshfunc);
    denseMSC.SetBuildArcGeometry(GInt::Vec3b(false, false, false));
    denseMSC.BuildArcs2DNoGeomDenseCache();
    const auto densePreHist = endpointHistogram(denseMSC);
    const INT_TYPE denseNodes = denseMSC.numNodes();
    const INT_TYPE denseArcs = denseMSC.numArcs();
    denseMSC.ComputeHierarchy(persLimit);

    const bool sameNodes = (genericNodes == denseNodes);
    const bool sameArcs = (genericArcs == denseArcs);
    const bool sameEndpoints = (genericPreHist == densePreHist);

    std::cout << "comparison: generic_nodes=" << genericNodes
              << " dense_nodes=" << denseNodes
              << " generic_arcs=" << genericArcs
              << " dense_arcs=" << denseArcs
              << " same_nodes=" << (sameNodes ? "true" : "false")
              << " same_arcs=" << (sameArcs ? "true" : "false")
              << " same_endpoint_histogram=" << (sameEndpoints ? "true" : "false")
              << std::endl;

    return 0;
}
