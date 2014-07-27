#include "myriaworld.h"
#include <iostream>
#include <fstream>
#include <limits>
#include <boost/geometry/strategies/transform.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/log/trivial.hpp>
#include "geo.hpp"
#include "SPHEmesh.hpp"

using myriaworld::country;
using myriaworld::polar2_point;
using myriaworld::polar2_polygon;
using myriaworld::polar2_multipolygon;

using myriaworld::cart2_point;
using myriaworld::cart2_polygon;
using myriaworld::cart2_multipolygon;

using myriaworld::cart3_point;
using myriaworld::cart3_polygon;
using myriaworld::cart3_multipolygon;

namespace bg = boost::geometry;

std::vector<country> read_countries(const std::string& filename){
        std::ifstream ifs(filename.c_str());
        unsigned n_poly, n_pts;
        double lat, lng;
        double minlat= std::numeric_limits<double>::infinity(),
               maxlat=-std::numeric_limits<double>::infinity();
        double minlng= std::numeric_limits<double>::infinity(),
               maxlng=-std::numeric_limits<double>::infinity();
        std::vector<country> countries;
        while(ifs){
            ifs >> n_poly;
            if(n_poly == 0)
                continue;
            country c;
            for(unsigned r=0; r< n_poly; r++){
                ifs >> n_pts;
                polar2_polygon poly;
                for(unsigned p=0; p<n_pts; p++){
                    ifs >> lat >> lng;
                    minlat = std::min(minlat, lat);
                    minlng = std::min(minlng, lng);
                    maxlat = std::max(maxlat, lat);
                    maxlng = std::max(maxlng, lng);
                    append(poly, polar2_point(lat, lng));
                }
                //BOOST_LOG_TRIVIAL(debug) << " - poly...";
                boost::geometry::correct(poly);
                //spherical_polygon poly2;
                //boost::geometry::simplify(poly, poly2, 0.001);
                //boost::geometry::correct(poly2);
                //if(poly2.outer().size() > 4)
                    //c.push_back(poly2);
                c.m_s2_polys.push_back(poly);
            }
            //assert(c.size() == n_poly);
            //BOOST_LOG_TRIVIAL(debug) << "# country...";
            if(c.m_s2_polys.size() > 0)
                countries.push_back(c);
            //assert(countries.back().size() == n_poly);
        }
        BOOST_LOG_TRIVIAL(info) << "countries -- minlat: "<< minlat << " maxlat: "<<maxlat;
        BOOST_LOG_TRIVIAL(info) << "countries -- minlng: "<< minlng << " maxlng: "<<maxlng;
        return countries;
}

std::vector<polar2_polygon> read_triangle_grid(const std::string& filename, int recdepth){
    BOOST_LOG_TRIVIAL(info) << "Reading triangles from `"<<filename<<"'...";
    std::vector<polar2_polygon>  trias;
    std::ifstream ifs(filename.c_str());
    double x0, y0, z0;
    double x1, y1, z1;
    double x2, y2, z2;
    double minlat= std::numeric_limits<double>::infinity(),
           maxlat=-std::numeric_limits<double>::infinity();
    double minlng= std::numeric_limits<double>::infinity(),
           maxlng=-std::numeric_limits<double>::infinity();
    while(ifs){
        ifs >> x0 >> y0 >> z0;
        ifs >> x1 >> y1 >> z1;
        ifs >> x2 >> y2 >> z2;
        if(!ifs)
            break;

        cart3_point p0(x0, y0, z0);
        cart3_point p1(x1, y1, z1);
        cart3_point p2(x2, y2, z2);
        polar2_point sp0, sp1, sp2;
        bg::strategy::transform::from_cartesian_3_to_spherical_equatorial_2<cart3_point, polar2_point> strategy;
        boost::geometry::transform(p0, sp0, strategy);
        boost::geometry::transform(p1, sp1, strategy);
        boost::geometry::transform(p2, sp2, strategy);
        double lat0 = sp0.get<0>()+2.1828, lng0 = sp0.get<1>();
        double lat1 = sp1.get<0>()+2.1828, lng1 = sp1.get<1>();
        double lat2 = sp2.get<0>()+2.1828, lng2 = sp2.get<1>();

        if(lat0 > 180) lat0 -= 360;
        if(lat1 > 180) lat1 -= 360;
        if(lat2 > 180) lat2 -= 360;

        //double eps = .01;
        //lat0 = lat0 - eps/2 + eps * drand48();
        //lat1 = lat1 - eps/2 + eps * drand48();
        //lat2 = lat2 - eps/2 + eps * drand48();

        //if(fabs(lat0) > 170) continue;
        //if(fabs(lat1) > 170) continue;
        //if(fabs(lat2) > 170) continue;

        minlat = std::min(minlat, lat0);
        minlng = std::min(minlng, lng0);
        maxlat = std::max(maxlat, lat0);
        maxlng = std::max(maxlng, lng0);

        polar2_polygon tria;
        append(tria, polar2_point(lat0, lng0));
        append(tria, polar2_point(lat1, lng1));
        append(tria, polar2_point(lat2, lng2));

        boost::geometry::correct(tria);
        if(bg::area(tria) > 0.000001)
            trias.emplace_back(tria);
    }
    BOOST_LOG_TRIVIAL(info) << "tria grid  -- minlat: "<< minlat << " maxlat: "<<maxlat;
    BOOST_LOG_TRIVIAL(info) << "tria grid  -- minlng: "<< minlng << " maxlng: "<<maxlng;
    return trias;
}

myriaworld::triangle_graph 
get_triangle_graph(int maxlevel){
    using namespace myriaworld;
    namespace bgi = boost::geometry::index;
    using boost::vertex_index_t;
    using boost::vertex_index;
    using boost::property_map;

    std::vector<cart3_polygon> ctriangles = triangulate_sphere(maxlevel);
    std::vector<polar2_polygon> triangles(ctriangles.size());

    for (unsigned int i = 0; i < triangles.size(); ++i)
    {
        boost::geometry::strategy::transform::from_cartesian_3_to_spherical_equatorial_2<cart3_point, polar2_point> strategy;
        boost::geometry::transform(ctriangles[i], triangles[i], strategy);
    }

    triangle_graph tg(triangles.size());

    boost::property_map<triangle_graph, vertex_index_t>::type index_map 
        = get(vertex_index, tg);

    boost::property_map<triangle_graph, vertex_pos_t>::type pos_map 
        = get(vertex_pos_t(), tg);

    boost::property_map<triangle_graph, n_shared_vert_t>::type n_shared_vert_map 
        = get(n_shared_vert_t(), tg);

    // the position, the vertex ID, the triangle ID
    typedef boost::tuple<cart2_point, unsigned, unsigned> value;
    bgi::rtree <value, bgi::quadratic<16> > rtree;
    std::vector<polar2_point> vertices;

    // add an edge between a and b if it is not yet in the set of edges
    auto add_edge_checked = [&](unsigned int tria, const polar2_point& a, const polar2_point& b){
        std::vector<value> res;
        using namespace myriaworld;

        cart2_point ca, cb;
        geo::cast_polar2_equatorial_2_cart_2(ca, a, 0.);
        geo::cast_polar2_equatorial_2_cart_2(cb, b, 0.);

        rtree.query(bgi::nearest(ca, 1), std::back_inserter(res));
        if(res.size() == 0 || geo::haversine_distance(vertices[res[0].get<1>()], a) > 0.001){
            vertices.push_back(a);
            rtree.insert({ca, vertices.size()-1, tria});
        }else{
            // we found another use of this vertex.
            int other = res.front().get<2>();
            auto e = boost::edge(other, tria, tg);
            if(!e.second){
                boost::add_edge(other, tria, tg);
                n_shared_vert_map[boost::edge(other, tria, tg).first]
                    = 1;
            }else{
                n_shared_vert_map[e.first]++;
            }
        }

        res.clear();
        rtree.query(bgi::nearest(cb, 1), std::back_inserter(res));
        if(res.size() == 0 || geo::haversine_distance(vertices[res[0].get<1>()], b) > 0.001){
            vertices.push_back(b);
            rtree.insert({cb, vertices.size()-1, tria});
        }else{
            // we found another use of this vertex.
            int other = res.front().get<2>();
            auto e = boost::edge(other, tria, tg);
            if(!e.second){
                boost::add_edge(other, tria, tg);
                n_shared_vert_map[boost::edge(other, tria, tg).first]
                    = 1;
            }else{
                n_shared_vert_map[e.first]++;
            }
        }
    };

    unsigned int tria_idx = 0;
    for(const auto& c : triangles){
        assert(c.outer().size() == 3);
        add_edge_checked(tria_idx, c.outer()[0], c.outer()[1]);
        add_edge_checked(tria_idx, c.outer()[1], c.outer()[2]);
        add_edge_checked(tria_idx, c.outer()[2], c.outer()[0]);
        ++tria_idx;
    }

    // remove all edges which have n_shared_vert != 2
    auto es = edges(tg);
    for(auto eit = es.first; eit!= es.second; eit++){
        if(n_shared_vert_map[*eit] < 2){
            boost::remove_edge(*eit, tg);
            continue;
        }
        assert(false);
    }
    return tg;
}
