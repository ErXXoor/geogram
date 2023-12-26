//
// Created by Hongbo on 12/25/23.
//

#ifndef GEOGRAM_LPCVTIS_H
#define GEOGRAM_LPCVTIS_H
#include <geogram/voronoi/integration_simplex.h>

namespace GEO{
    class GEOGRAM_API LpCVTIS : public IntegrationSimplex {
    public:
    LpCVTIS(const Mesh &mesh, bool volumetric, index_t nb_frames,
            index_t nb_comp_per_frame, const double *frames);

    ~LpCVTIS() override = default;

    const double* vertex_ptr(index_t v) const;

    void set_dim(unsigned int dim) { m_dim = dim; }
    double eval(
            index_t center_vertex_index,
            const GEOGen::Vertex& v0,
            const GEOGen::Vertex& v1,
            const GEOGen::Vertex& v2,
            index_t t,
            index_t t_adj = index_t(-1),
            index_t v_adj = index_t(-1)
            ) override;

    private:
    unsigned int m_dim;

    };

}
#endif //GEOGRAM_LPCVTIS_H
