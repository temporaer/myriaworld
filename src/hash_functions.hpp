#ifndef __HASH_FUNCTIONS_HPP_876525__
#     define __HASH_FUNCTIONS_HPP_876525__

#include "myriaworld.h"
namespace myriaworld{
    const triangle_graph* g_tg;
    std::size_t hash_value(const 
            boost::graph_traits<triangle_graph>::edge_iterator::value_type&
            e){
        std::size_t seed=0;
        boost::property_map<triangle_graph, boost::edge_weight_t>::const_type weightmap = boost::get(boost::edge_weight, *g_tg);
        boost::hash_combine(seed, weightmap[e]);
        return seed;
    }
    std::size_t hash_value(const country_bit& c){
        std::size_t seed=0;
        boost::hash_combine(seed, c.m_scalerank);
        boost::hash_combine(seed, c.m_mapcolor);
        boost::hash_combine(seed, c.m_s2_poly.outer());
        return seed;
    }
    std::size_t hash_value(const triavertex_descriptor& v){
        std::size_t seed=0;
        boost::property_map<triangle_graph, vertex_pos_t>::const_type pos_map 
            = get(vertex_pos_t(), *g_tg);
        boost::property_map<triangle_graph, vertex_fracfilled_t>::const_type frac_filled 
            = get(vertex_fracfilled_t(), *g_tg);
        boost::property_map<triangle_graph, country_bits_t>::const_type bs_map 
            = get(country_bits_t(), *g_tg);
        boost::hash_combine(seed,
                boost::hash_range(pos_map[v].m_s2_poly.outer().begin(),
                                  pos_map[v].m_s2_poly.outer().end()));
        boost::hash_combine(seed, frac_filled[v]);
        boost::hash_combine(seed, frac_filled[v]);
        boost::hash_combine(seed,
                boost::hash_range(bs_map[v].begin(),
                                  bs_map[v].end()));
                             
        return seed;
    }

    std::size_t hash_value(const triangle_graph& g){
        std::size_t seed = 0;
        auto vs = vertices(g);
        auto es = edges(g);
        g_tg = &g;
        boost::property_map<triangle_graph, boost::edge_weight_t>::const_type weightmap = boost::get(boost::edge_weight, *g_tg);
        boost::hash_combine(seed,
                boost::hash_range(vs.first, vs.second));
        for(; es.first != es.second; es.first++)
            boost::hash_combine(seed, weightmap[*es.first]);
        return seed;
    }
    std::size_t hash_value(
            const std::function<void(myriaworld::triangle_graph&, const std::vector<myriaworld::country>&)>&){
        return 0;
    }
    std::size_t hash_value(
            const std::function<void(myriaworld::triangle_graph&)>&){
        return 0;
    }
}

#endif /* __HASH_FUNCTIONS_HPP_876525__ */
