#ifndef __MYRIAWORLD_HPP__
#     define __MYRIAWORLD_HPP__
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace myriaworld{
    namespace bg = boost::geometry;
    typedef bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree> > polar2_point;
    typedef bg::model::point<double, 2, bg::cs::cartesian> cart2_point;
    typedef bg::model::point<double, 3, bg::cs::cartesian> cart3_point;

    typedef bg::model::polygon<cart2_point, true, false> cart2_polygon;
    typedef bg::model::polygon<cart3_point, true, false> cart3_polygon;
    typedef bg::model::polygon<polar2_point, true, false> polar2_polygon;

    typedef bg::model::multi_polygon<polar2_polygon> polar2_multipolygon;
    typedef bg::model::multi_polygon<cart3_polygon> cart3_multipolygon;
    typedef bg::model::multi_polygon<cart2_polygon> cart2_multipolygon;


    typedef std::pair<int, int> edge_t;

    struct city{
        std::string m_name;
        int m_rank;
        cart3_point m_c3_location;   // 3d coordinate on unit sphere
        polar2_point m_s2_location;  // coordinate on map
        cart2_point m_c2_location;   // as read from shapefile
    };

    struct country{
        std::string m_name;          // the name of this country
        polar2_multipolygon m_s2_polys;  // all polygons which constitute the country
    };

    struct triangle{
        cart3_polygon m_c3_poly;   // 3d coordinate on unit sphere
        cart2_polygon m_c2_poly;   // coordinate on map
        polar2_polygon m_s2_poly;  // as read produced by sphere triangulation
    };

    struct country_bit{
        boost::shared_ptr<country> m_country;
        cart3_polygon m_c3_poly;   // 3d coordinate on unit sphere
        cart2_polygon m_c2_poly;   // coordinate on map
        polar2_polygon m_s2_poly;  // as read from shapefile
        country_bit(){}
        country_bit(const polar2_polygon& sp):m_s2_poly(sp){}
    };

    struct shared_edge_property{
        // the "lower" vertex is the one with the lower index in the graph
        
        // the "shared" vars are indices referring to the ones which are in
        // source and target vertex (=triangle)
        unsigned char shared_l0, shared_l1, shared_h0, shared_h1;

        // the vertex which is not present in both triangles
        // (index in triangle-polygon)
        unsigned char single_l, single_h;

        template<class T>
        inline shared_edge_property flipped_if_needed(const T& a, const T& b){
            if(&a < &b)
                return *this;
            shared_edge_property ret(*this);
            std::swap(ret.shared_l0, ret.shared_h0);
            std::swap(ret.shared_l1, ret.shared_h1);
            std::swap(ret.single_l, ret.single_h);
            return ret;
        }

        bool is_cut;
    };

    typedef std::vector<country_bit> country_bit_vec;
    typedef std::vector<city> city_vec;

    struct vertex_pos_t        { typedef boost::vertex_property_tag kind; };
    struct vertex_area_t       { typedef boost::vertex_property_tag kind; };
    struct vertex_fracfilled_t { typedef boost::vertex_property_tag kind; };
    struct country_bits_t      { typedef boost::vertex_property_tag kind; };
    struct cities_t            { typedef boost::vertex_property_tag kind; };
    struct n_shared_vert_t     { typedef boost::edge_property_tag kind;   };
    struct shared_edge_t       { typedef boost::edge_property_tag kind;   };

    typedef boost::property<boost::edge_weight_t, double,
            boost::property<n_shared_vert_t, int,
            boost::property<shared_edge_t, shared_edge_property > > >
        TriangleEdgeProperty;
    typedef boost::property<vertex_pos_t, triangle,
            boost::property<vertex_area_t, double,
            boost::property<vertex_fracfilled_t, double,
            boost::property<country_bits_t, country_bit_vec,
            boost::property<cities_t, city_vec> > > > >
                TriangleVertexProperty;

    typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS,
            TriangleVertexProperty, TriangleEdgeProperty> triangle_graph;

    typedef boost::graph_traits<triangle_graph>::vertex_descriptor triavertex_descriptor;
    typedef boost::graph_traits<triangle_graph>::edge_descriptor triaedge_descriptor;
}
#endif /* __MYRIAWORLD_HPP__ */
