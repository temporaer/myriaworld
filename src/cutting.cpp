#include "cutting.hpp"
#include "geo.hpp"
#include <boost/graph/prim_minimum_spanning_tree.hpp>

namespace myriaworld
{
    void determine_cuttings(triangle_graph& g){
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
    }

    void flatten(triangle_graph& g, double fact){
        std::vector < boost::graph_traits<triangle_graph>::vertex_descriptor >
            p(boost::num_vertices(g));
        boost::property_map<triangle_graph, vertex_pos_t>::type pos_map 
            = get(vertex_pos_t(), g);
        boost::property_map<triangle_graph, shared_edge_t>::type se_map 
            = get(shared_edge_t(), g);


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
            open_list.push_back(0);

            // first 3d point on map
            pos_map[0].m_mappos = pos_map[0].m_c3_poly;


            while(open_list.size() > 0){
                // current triangle in globe coordinates
                int cur_tria = open_list.back();
                open_list.pop_back();
                printf("flattening triangle %d\n", cur_tria);

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

                    open_list.push_back(targetv);
                }
            }
        }
    }
}
