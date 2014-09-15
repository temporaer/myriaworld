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

    // true == clockwise, closed
    // Google S2 requires CCW!
    // S2 polygons consist of loops, and loops are /implicitly/ closed.
    typedef bg::model::polygon<cart2_point  , false, false> cart2_polygon;
    typedef bg::model::polygon<cart3_point  , false, false> cart3_polygon;
    typedef bg::model::polygon<polar2_point , false, false> polar2_polygon;

    typedef bg::model::multi_polygon<polar2_polygon> polar2_multipolygon;
    typedef bg::model::multi_polygon<cart3_polygon> cart3_multipolygon;
    typedef bg::model::multi_polygon<cart2_polygon> cart2_multipolygon;


    typedef std::pair<int, int> edge_t;

    struct city{
        std::string m_name;
        int m_rank;
        cart3_point m_c3_location;   // 3d coordinate on unit sphere
        polar2_point m_s2_location;  // as read from shapefile
        cart3_point m_mappos;   // as read from shapefile
        template<class Archive>
        void serialize(Archive& ar, const unsigned long ){
            ar & m_name & m_rank & m_c3_location & m_s2_location & m_mappos;
        }
    };

    struct country{
        std::string m_name;          // the name of this country
        double m_mapcolor;
        double m_scalerank;
        polar2_multipolygon m_s2_polys;  // all polygons which constitute the country
        template<class Archive>
        void serialize(Archive& ar, const unsigned long ){
            ar & m_name & m_mapcolor & m_scalerank & m_s2_polys;
        }
    };

    struct triangle{
        cart3_polygon m_c3_poly;   // 3d coordinate on unit sphere
        cart3_polygon m_mappos;   // coordinate on map
        polar2_polygon m_s2_poly;  // as read produced by sphere triangulation
        template<class Archive>
        void serialize(Archive& ar, const unsigned long ){
            ar & m_c3_poly & m_mappos & m_s2_poly;
        }
    };

    struct country_bit{
        cart3_polygon m_c3_poly;   // 3d coordinate on unit sphere
        double m_mapcolor;
        double m_scalerank;
        std::string m_name;
        cart3_polygon m_mappos;   // coordinate on map
        polar2_polygon m_s2_poly;  // as read from shapefile
        country_bit(const polar2_polygon& sp):m_s2_poly(sp){}
        template<class Archive>
        void serialize(Archive& ar, const unsigned long ){
            ar & m_c3_poly & m_mappos & m_s2_poly & m_mappos & m_scalerank & m_mapcolor & m_scalerank & m_name;
        }
        country_bit(){}
        country_bit(const country& c)
            : m_mapcolor(c.m_mapcolor),
              m_scalerank(c.m_scalerank),
              m_name(c.m_name)
        {}
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
            if(a < b)
                return *this;
            shared_edge_property ret(*this);
            std::swap(ret.shared_l0, ret.shared_h0);
            std::swap(ret.shared_l1, ret.shared_h1);
            std::swap(ret.single_l, ret.single_h);
            return ret;
        }

        bool is_cut;
        template<class Archive>
        void serialize(Archive& ar, const unsigned long ){
            ar & shared_l0 & shared_l1 & shared_h0 & shared_h1;
            ar & single_l & single_h;
            ar & is_cut;
        }
    };

    typedef std::vector<country_bit> country_bit_vec;
    typedef std::vector<city> city_vec;

    struct vertex_pos_t        { typedef boost::vertex_property_tag kind; };
    struct vertex_area_t       { typedef boost::vertex_property_tag kind; };
    struct vertex_centroid_t   { typedef boost::vertex_property_tag kind; };
    struct vertex_fracfilled_t { typedef boost::vertex_property_tag kind; };
    struct country_bits_t      { typedef boost::vertex_property_tag kind; };
    struct cities_t            { typedef boost::vertex_property_tag kind; };
    struct parent_t            { typedef boost::vertex_property_tag kind; };
    struct n_shared_vert_t     { typedef boost::edge_property_tag kind;   };
    struct shared_edge_t       { typedef boost::edge_property_tag kind;   };

    typedef boost::property<boost::edge_weight_t, double,
            boost::property<n_shared_vert_t, int,
            boost::property<shared_edge_t, shared_edge_property > > >
        TriangleEdgeProperty;
    typedef boost::property<vertex_pos_t, triangle,
            boost::property<vertex_area_t, double,
            boost::property<vertex_centroid_t, polar2_point,
            boost::property<vertex_fracfilled_t, double,
            boost::property<country_bits_t, country_bit_vec,
            boost::property<parent_t, size_t,
            boost::property<cities_t, city_vec> > > > > > >
                TriangleVertexProperty;

    typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS,
            TriangleVertexProperty, TriangleEdgeProperty> triangle_graph;

    typedef boost::graph_traits<triangle_graph>::vertex_descriptor triavertex_descriptor;
    typedef boost::graph_traits<triangle_graph>::edge_descriptor triaedge_descriptor;
}
#endif /* __MYRIAWORLD_HPP__ */
