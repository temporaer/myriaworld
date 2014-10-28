#include "cutting.hpp"
#include "geo.hpp"
#include <boost/log/trivial.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace myriaworld
{
    triangle_graph determine_cuttings(triangle_graph& g){
        std::vector < boost::graph_traits<triangle_graph>::vertex_descriptor >
            p(boost::num_vertices(g));
        boost::property_map<triangle_graph, vertex_pos_t>::type pos_map 
            = get(vertex_pos_t(), g);
        boost::property_map<triangle_graph, parent_t>::type pred_map 
            = get(parent_t(), g);
        boost::property_map<triangle_graph, shared_edge_t>::type se_map 
            = get(shared_edge_t(), g);

        // 4a. find a vertex at the "back" of the globe (somewhere around samoa ^^)
        triavertex_descriptor max_vertex;
        {
            double max_vertex_val = 1e12;
            auto vs = vertices(g);
            for(auto vit = vs.first; vit!=vs.second; vit++){
                auto pos = geo::spherical_centroid(pos_map[*vit].m_s2_poly);
                double dist = geo::haversine_distance(pos, polar2_point(88, 46));
                if(dist < max_vertex_val){
                    max_vertex_val = dist;
                    max_vertex = *vit;
                }
            }
        }
        boost::prim_minimum_spanning_tree(g, &p[0], boost::root_vertex(max_vertex));
        auto es          = boost::edges(g);
        pred_map[0] = 0;
        for(auto eit = es.first; eit != es.second; eit++){
            const auto sourcev = source(*eit, g);
            const auto targetv = target(*eit, g);
            bool is_part_of_mst = false;
            if(p[sourcev] == targetv){
                is_part_of_mst = true;
                pred_map[sourcev] = targetv;
            }
            if(p[targetv] == sourcev){
                is_part_of_mst = true;
                pred_map[targetv] = sourcev;
            }
            se_map[*eit].is_cut = !is_part_of_mst;
        }
        return g;
    }

    /// return same triangle, but in the xy plane.
    cart3_polygon tria3dto2d(const cart3_polygon& p){
        using namespace geo;
        cart3_polygon q;

        cart3_point d0 = difference(p.outer()[1], p.outer()[0]);
        cart3_point d1 = difference(p.outer()[2], p.outer()[0]);
        //cart3_point d2 = difference(p.outer()[2], p.outer()[1]);
        //cart3_point n0 = cross(d0, d1);

        double l0 = norm(d0);
        double l1 = norm(d1);
        //double l2 = norm(d2);
        double costheta = fabs(dot(divide(d0,norm(d0)), divide(d1,norm(d1))));

        q.outer().push_back(cart3_point(0, 0, 0));
        q.outer().push_back(cart3_point(l0, 0, 0));
        q.outer().push_back(cart3_point(costheta*l0, sin(acos(costheta))*l1, 0));
        return q;
    }

    void flatten(triangle_graph& g, double fact, double clat, double clon){
        std::vector < boost::graph_traits<triangle_graph>::vertex_descriptor >
            p(boost::num_vertices(g));
        boost::property_map<triangle_graph, vertex_pos_t>::type pos_map 
            = get(vertex_pos_t(), g);
        boost::property_map<triangle_graph, shared_edge_t>::type se_map 
            = get(shared_edge_t(), g);
        boost::property_map<triangle_graph, country_bits_t>::type cb_map 
            = get(country_bits_t(), g);
        boost::property_map<triangle_graph, vertex_centroid_t>::type centroid_map
            = get(vertex_centroid_t(), g);


        {
            std::vector<triavertex_descriptor> closed_list;

            // transform all spherical coordinates to cart3 coordinates
            auto vs = boost::vertices(g);
            for(auto vit = vs.first; vit != vs.second; vit++){
                using namespace boost::geometry::strategy::transform;
                from_spherical_equatorial_2_to_cartesian_3<polar2_point, cart3_point> sph2cart3;
                bg::transform(pos_map[*vit].m_s2_poly, pos_map[*vit].m_c3_poly, sph2cart3);
            }

            using namespace boost;
            std::vector<triavertex_descriptor> open_list;

            int start_idx = 0;
            {
                int s = num_vertices(g);
                double min_d = 1E9;
                polar2_point cmp(clon,clat);
                BOOST_LOG_TRIVIAL(warning) << "   clon " << clon << " clat " <<clat;
                for(unsigned int i=0; i< s; i++){
                    double dist = geo::haversine_distance(centroid_map[i], cmp);
                    
                    if(dist < min_d) {
                        min_d = dist;
                        start_idx = i;
                    }
                    //BOOST_LOG_TRIVIAL(warning) << "   c "<<  centroid_map[i].get<0>() << ", " << centroid_map[i].get<1>();
                    //BOOST_LOG_TRIVIAL(warning) << "   i " << i << " si " << start_idx << " d " << min_d;
                }
            }
            BOOST_LOG_TRIVIAL(warning) << "start_idx " << start_idx;

            open_list.push_back(start_idx);

            // first 3d point on map
            //pos_map[start_idx].m_mappos = tria3dto2d(pos_map[start_idx].m_c3_poly);
            {
                boost::geometry::transform(pos_map[start_idx].m_c3_poly, pos_map[start_idx].m_s2_poly);
                boost::geometry::transform(geo::rotated(pos_map[start_idx].m_s2_poly, 
                            -centroid_map[start_idx].get<0>(),
                            -centroid_map[start_idx].get<1>()), 
                        pos_map[start_idx].m_mappos);
            }


            while(open_list.size() > 0){
                // current triangle in globe coordinates
                int cur_tria = open_list.back();
                open_list.pop_back();

                closed_list.push_back(cur_tria);
                //if(closed_list.size() > 10)
                    //break;

                // find neighboring triangles
                auto es          = boost::out_edges(cur_tria, g);
                for(auto eit = es.first; eit != es.second; eit++){
                    const auto sourcev = source(*eit, g);
                    const auto targetv = target(*eit, g);

                    // avoid doing any edge twice
                    if(std::find(closed_list.begin(), closed_list.end(), targetv) != closed_list.end())
                        continue;

                    if(se_map[*eit].is_cut)
                        continue;

                    // calculate angle between them on original globe
                    // put new point onto same plane as current triangle
                    pos_map[targetv].m_mappos 
                        = geo::flatten(pos_map[sourcev].m_c3_poly, pos_map[targetv].m_c3_poly, 
                            pos_map[sourcev].m_mappos, se_map[*eit].flipped_if_needed(sourcev, targetv), fact);
                    
                    for(auto& b : cb_map[targetv]){
                        bg::strategy::transform::from_spherical_equatorial_2_to_cartesian_3<polar2_point, cart3_point> strategy;
                        //bg::correct(b.m_s2_poly);
                        boost::geometry::transform(b.m_s2_poly, b.m_c3_poly, strategy);
                        //bg::correct(b.m_c3_poly);
                        b.m_mappos = geo::tria2tria(
                                pos_map[targetv].m_c3_poly,
                                pos_map[targetv].m_mappos,
                                b.m_c3_poly);
                    }

                    open_list.push_back(targetv);
                }
            }
        }
        if(0){
            // merging close-by points into one
            namespace bgi = boost::geometry::index;
            typedef boost::tuple<cart3_point, cart3_point*, unsigned long> value;
            bgi::rtree <value, bgi::quadratic<16> > rtree;
            std::vector<bool> adjusted;
            adjusted.reserve(1000000);
            auto vs = boost::vertices(g);
            unsigned long cnt = 0;
            for(auto vit = vs.first; vit != vs.second; vit++){
                for(auto& bits : cb_map[*vit]){
                    for(auto& pos : bits.m_mappos.outer()){
                        rtree.insert({pos, &pos, cnt++});
                        adjusted.push_back(false);
                    }
                }
            }
            cnt = 0;
            std::vector<value> qres;
            std::vector<cart3_point*> actual_neighbors;
            actual_neighbors.reserve(1000); // plenty
            qres.reserve(5000); // plenty
            for(auto vit = vs.first; vit != vs.second; vit++){
                for(auto& bits : cb_map[*vit]){
                    for(auto& pos : bits.m_mappos.outer()){
                        if(adjusted[cnt++])
                            continue;
                        const static double eps = 0.003;
                        auto box = boost::geometry::model::box<cart3_point>(
                                cart3_point(pos.get<0>() - eps/2., pos.get<1>() - eps/2., pos.get<2>() - eps/2.),
                                cart3_point(pos.get<0>() + eps/2., pos.get<1>() + eps/2., pos.get<2>() + eps/2.)
                                );
                        qres.clear();
                        rtree.query(bgi::intersects(box), std::back_inserter(qres));
                        cart3_point sum {0, 0, 0};
                        int wsum = 0;
                        actual_neighbors.clear();
                        for(auto& qp : qres){
                            if(boost::geometry::distance(qp.get<0>(), pos) < eps){
                                sum = geo::sum(qp.get<0>(), sum);
                                wsum ++;
                                actual_neighbors.push_back(qp.get<1>());
                                adjusted[qp.get<2>()] = true;
                            }
                        }
                        sum = geo::divide(sum, (double) wsum);
                        for(auto& qp : actual_neighbors){
                            *qp = sum;
                        }
                    }
                }
            }
        }
    }
}
