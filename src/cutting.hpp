#ifndef __CUTTING_HPP__
#     define __CUTTING_HPP__
#include "myriaworld.h"

namespace myriaworld
{
    triangle_graph determine_cuttings(triangle_graph& g);
    void flatten(triangle_graph& g, double fact, double clat, double clon, double roll);
}
#endif /* __CUTTING_HPP__ */
