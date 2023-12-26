//
// Created by Hongbo on 12/25/23.
//
#include "LpCVTIS.h"
#include <geogram/voronoi/generic_RVD.h>
#include <geogram/basic/smart_pointer.h>
namespace GEO {
    LpCVTIS::LpCVTIS(const Mesh &mesh,
                      bool volumetric,
                      index_t nb_frames,
                      index_t nb_comp_per_frame,
                      const double *frames) :
                      IntegrationSimplex(mesh,
                                        volumetric,
                                        nb_frames,
                                        nb_comp_per_frame,
                                        frames){
    }

    const double* LpCVTIS::vertex_ptr(index_t i) const {
        return points_ + i * points_stride_;
    }

    double LpCVTIS::eval(
            index_t v,
            const GEOGen::Vertex &p1,
            const GEOGen::Vertex &p2,
            const GEOGen::Vertex &p3,
            index_t t,
            index_t t_adj,
            index_t v_adj
    ) {
        auto p0 = vertex_ptr(v);
        double t_area = Geom::triangle_area<double>(p1, p2, p3, m_dim);
        double cur_f = 0.0;
        for(index_t c = 0; c < m_dim; c++) {
            double u0 = p0[c] - p1[c];
            double u1 = p0[c] - p2[c];
            double u2 = p0[c] - p3[c];
            cur_f += u0 * u0;
            cur_f += u1 * (u0 + u1);
            cur_f += u2 * (u0 + u1 + u2);
        }
        double f = t_area * cur_f / 6.0;

        for(index_t c = 0; c < m_dim; c++) {
            double Gc = (1.0 / 3.0) * (p1[c] + p2[c] + p3[c]);
            g_[m_dim * v + c] += (2.0 * t_area) * (p0[c] - Gc);
        }
        return f;
    }
}

