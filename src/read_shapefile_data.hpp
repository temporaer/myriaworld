#ifndef __READ_SHAPEFILE_DATA_HPP__
#     define __READ_SHAPEFILE_DATA_HPP__
#include "myriaworld.h"

namespace myriaworld
{

    std::vector<myriaworld::country> read_countries(const std::string& filename);

    std::vector<polar2_polygon> read_triangle_grid(const std::string& filename);
    myriaworld::triangle_graph get_triangle_graph(int maxlevel);

    myriaworld::triangle_graph smooth_landmass(myriaworld::triangle_graph& g, double sigma);
    myriaworld::triangle_graph determine_edge_weights(myriaworld::triangle_graph& g, double wlat, double wlon);

}
#endif /* __READ_SHAPEFILE_DATA_HPP__ */
