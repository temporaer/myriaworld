#ifndef __READ_SHAPEFILE_DATA_HPP__
#     define __READ_SHAPEFILE_DATA_HPP__
#include "myriaworld.h"

std::vector<myriaworld::country> read_countries(const std::string& filename);

myriaworld::triangle_graph read_triangles(const std::string& filename);

#endif /* __READ_SHAPEFILE_DATA_HPP__ */
