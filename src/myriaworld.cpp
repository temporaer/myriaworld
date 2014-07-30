#include <fstream>
#include "myriaworld.h"
#include "read_shapefile_data.hpp"
#include "cutting.hpp"
#include "country2tria.hpp"
#include "geo.hpp"
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/serialization/vector.hpp>
//#include <boost/geometry/index/detail/serialization.hpp>
#include "boost_geo_serialize.hpp"
#include "cached_function.hpp"

void write_s2centroids(std::string filename, myriaworld::triangle_graph& g){
    BOOST_LOG_TRIVIAL(info) << "Writing to file: " << filename;
    using myriaworld::polar2_point;
    using myriaworld::cart3_point;
    using myriaworld::triangle_graph;
    using myriaworld::vertex_pos_t;
    using myriaworld::shared_edge_t;
    using boost::edge_weight_t;
    using boost::edge_weight;
    std::vector<polar2_point> centroids(boost::num_vertices(g));
    boost::property_map<triangle_graph, vertex_pos_t>::type pos_map = get(vertex_pos_t(), g);
    boost::property_map<triangle_graph, edge_weight_t>::type weight_map = get(edge_weight, g);
    boost::property_map<triangle_graph, shared_edge_t>::type se_map = get(shared_edge_t(), g);
    auto vs = boost::vertices(g);
    for(auto vit = vs.first; vit != vs.second; vit++){
        centroids[*vit] = myriaworld::geo::spherical_centroid(pos_map[*vit].m_s2_poly);
    }

    std::ofstream ofs(filename.c_str());
    auto optional_print_s2 = [&](double color, const polar2_point& a, const polar2_point& b){
        if(myriaworld::geo::crosses_dateline(a, b))
            return;
        ofs << color << " " << a.get<0>() << " " << a.get<1>()
            << " "          << b.get<0>() << " " << b.get<1>() << std::endl;
    };
    auto optional_print_c2 = [&](double color, const cart3_point& a, const cart3_point& b){
        ofs << color << " " << a.get<0>() << " " << a.get<1>()
            << " "          << b.get<0>() << " " << b.get<1>() << std::endl;
    };

    //for(auto vit = vs.first; vit != vs.second; vit++){
    //    //auto t = pos_map[*vit].m_mappos;
    //    //if(t.outer().size() != 3)
    //        //continue;
    //    //optional_print_c2(1, t.outer()[0], t.outer()[1]);
    //    //optional_print_c2(1, t.outer()[0], t.outer()[2]);
    //    //optional_print_c2(1, t.outer()[1], t.outer()[2]);
    //    auto t = pos_map[*vit].m_s2_poly;
    //    optional_print_s2(1, t.outer()[0], t.outer()[1]);
    //    optional_print_s2(1, t.outer()[0], t.outer()[2]);
    //    optional_print_s2(1, t.outer()[1], t.outer()[2]);
    //}
    auto es = boost::edges(g);
    std::cout << "n edges: " << boost::num_edges(g) << std::endl;
    std::cout << "n vert: " << boost::num_vertices(g) << std::endl;
    for(auto eit = es.first; eit != es.second; eit++){
        auto sourcev = boost::source(*eit, g);
        auto targetv = boost::target(*eit, g);
        auto se = myriaworld::geo::determine_shared_vertices(pos_map[sourcev].m_s2_poly, pos_map[targetv].m_s2_poly);
        //auto s = centroids[boost::source(*eit, g)];
        //auto t = centroids[boost::target(*eit, g)];
        auto p = pos_map[sourcev];
        //if(se_map[*eit].is_cut)
            //optional_print_s2(1 + weight_map[*eit], s, t);
        //if(p.m_mappos.outer().size() != 3){
        //    BOOST_LOG_TRIVIAL(info) << "wrong number of vertices: " << p.m_mappos.outer().size();
        //    continue;
        //}
        //auto s = p.m_s2_poly.outer()[se.shared_l0];
        //auto t = p.m_s2_poly.outer()[se.shared_l1];
        auto s = p.m_mappos.outer()[se.shared_l0];
        auto t = p.m_mappos.outer()[se.shared_l1];
        if(!se_map[*eit].is_cut)
            ;//optional_print_c2(weight_map[*eit], s, t);
        else
            optional_print_c2(weight_map[*eit], s, t);
    }
}

//template<class Func, class T, class... Args>
//inline 
//T& side_effect_wrap(Func f, T& t, Args&&... args){
    //f(t, args...);
    //return t;
//}

myriaworld::triangle_graph side_effect_wrap(
        std::function<void(myriaworld::triangle_graph&, const std::vector<myriaworld::country>&)> func,
        myriaworld::triangle_graph& g, const std::vector<myriaworld::country>& c){
    func(g, c);
    return g;
}


int main(){
    using namespace myriaworld;
    auto countries = read_countries("../../myriaworld/easy.txt");
    auto g = get_triangle_graph(5);
    fscache::function_cache c;
    c("country2tria", 0, side_effect_wrap, country2tria, g, countries);
    determine_edge_weights(g);
    determine_cuttings(g);
    write_s2centroids("triagrid.txt", g);
    flatten(g, 0.);
}
