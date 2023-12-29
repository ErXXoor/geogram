//
// Created by Hongbo on 12/25/23.
//

#ifndef GEOGRAM_LPCVTIS_H
#define GEOGRAM_LPCVTIS_H
#include <geogram/voronoi/integration_simplex.h>

namespace GEO{
    class GEOGRAM_API LpCVTIS : public IntegrationSimplex {
    public:
    LpCVTIS(const Mesh &mesh, bool volumetric, unsigned int dim, unsigned int degree);

    ~LpCVTIS() override = default;

    const double* vertex_ptr(index_t v) const;
    double eval(
            index_t center_vertex_index,

            const GEOGen::Vertex& v0,
            const GEOGen::Vertex& v1,
            const GEOGen::Vertex& v2,
            index_t t,
            index_t t_adj = index_t(-1),
            index_t v_adj = index_t(-1)
            ) override;

    double grad_tri(const vec3& U1, const vec3& U2, const vec3& U3,
                        vec3& dTdU1, vec3& dTdU2, vec3& dTdU3);
    //Utils
    void vecmul(const double* p1, const double* p2, double* to);
    void vecmul(const double* p1, const double* p2, const double* p3, double* to);
    double vecbar(const double* p1);
    void vecmadd(double s, const double* p1, const double* p2, const double* p3, double* to);
    void vecmadd(double s, const vec3& p1, double t, const vec3& p2, vec3& to);
    void matTvecmul(const mat3& M, const vec3& U, vec3& V);

    private:
    unsigned int m_dim;
    unsigned int m_degree;
    unsigned int nb_coeffs;
    unsigned int nb_dcoeffs;
    std::vector<std::vector<unsigned int>> E_pow;
    std::vector<std::vector<std::vector<unsigned int>>> dE_pow;

    };

}
#endif //GEOGRAM_LPCVTIS_H
