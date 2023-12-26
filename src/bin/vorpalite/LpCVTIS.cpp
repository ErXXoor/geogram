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
    double LpCVTIS::eval(
            index_t v,
            const GEOGen::Vertex &v0,
            const GEOGen::Vertex &v1,
            const GEOGen::Vertex &v2,
            index_t t,
            index_t t_adj,
            index_t v_adj
    ) {
         auto v3 = mesh_.vertices.point_ptr(v);
         double t_area = Geom::triangle_area(*v0.point_ptr(), *v1.point_ptr(), *v2.point_ptr());
         return 0;
    }
}

