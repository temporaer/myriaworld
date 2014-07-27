#include <boost/geometry/index/rtree.hpp>
#include <boost/log/trivial.hpp>
#include "myriaworld.h"
#include "country2tria.hpp"
#include "geo.hpp"

namespace myriaworld
{
    void country2tria(triangle_graph& g, const std::vector<country>& countries){
        namespace bgi = boost::geometry::index;

        std::vector<polar2_polygon> flat_countries;

        for(const auto& c : countries)
            for(const auto& poly : c.m_s2_polys)
                flat_countries.push_back(poly);

        std::vector<bool> done_checks(boost::num_vertices(g) * flat_countries.size());

        typedef bg::model::box<cart2_point> box;
        typedef std::pair<box, unsigned> value;

        unsigned int n_intersection_checks = 0;
        unsigned int box_lat = 145, box_lng = 70;
        for(int lat_off = -40; lat_off <= 40; lat_off += 20){
            for(int lng_off = -40; lng_off <= 40; lng_off += 20){
                bgi::rtree <value, bgi::quadratic<16> > country_idx;
                {   unsigned idx=0;
                    box b;
                    for(const auto& poly2 : flat_countries){
                        polar2_polygon poly = geo::rotated(poly2, lat_off, lng_off);
                        if(geo::inside_box(poly, box_lat, box_lng)){
                            cart2_polygon cp;
                            geo::cast_polar2_equatorial_2_cart_2(cp, poly, 0);
                            bg::correct(cp);
                            bg::envelope(cp, b);
                            country_idx.insert(std::make_pair(b, idx));
                        }
                        idx++;
                    }
                }

                unsigned int tria_id = 0;
                auto verts = boost::vertices(g);
                boost::property_map<triangle_graph, vertex_pos_t>::type pos_map
                    = get(vertex_pos_t(), g);
                boost::property_map<triangle_graph, country_bits_t>::type bits_map
                    = get(country_bits_t(), g);
                for(auto vit = verts.first; vit != verts.second; vit++){
                    polar2_polygon tria2 = pos_map[*vit].m_s2_poly;
                    polar2_polygon tria = geo::rotated(tria2, lat_off, lng_off);
                    if( !geo::inside_box(tria, box_lat, box_lng)){
                        tria_id++;
                        continue;
                    }

                    cart2_polygon etria;
                    geo::cast_polar2_equatorial_2_cart_2(etria, tria, 0);
                    bg::correct(etria);

                    std::vector<value> result_c;
                    country_idx.query(bgi::intersects(etria),
                            std::back_inserter(result_c));
                    for(const auto& cres : result_c){
                        if(done_checks[tria_id * flat_countries.size() + cres.second])
                            continue;
                        done_checks[tria_id * flat_countries.size() + cres.second] = true;

                        const polar2_polygon& c2 = flat_countries[cres.second];
                        polar2_polygon c = geo::rotated(c2, lat_off, lng_off);

                        cart2_polygon ep;
                        geo::cast_polar2_equatorial_2_cart_2(ep, c, 0.);
                        boost::geometry::correct(ep);

                        std::vector<cart2_polygon> output;
                        try{
                            bg::intersection(etria, ep, output);
                            n_intersection_checks ++;
                        }catch(...){
                            BOOST_LOG_TRIVIAL(warning) << "Error in boost intersection";
                            continue;
                        }
                        for(auto& rpoly : output)
                        {
                            boost::geometry::correct(rpoly); // necessary bc of BOOST_GEOMETRY_OVERLAY_NO_THROW
                            polar2_polygon sp;
                            geo::cast_cart_2_polar2_equatorial_2(sp, rpoly, 0.);
                            boost::geometry::correct(sp);
                            bits_map[*vit].push_back(country_bit(geo::rotated(sp, -lat_off, -lng_off)));
                        }
                    }
                    ++tria_id;
                }
            }
        }
    }

}
