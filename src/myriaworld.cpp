#include <fstream>
#include "myriaworld.h"
#include "read_shapefile_data.hpp"
#include "cutting.hpp"
#include "country2tria.hpp"
#include "geo.hpp"
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
//#include <boost/geometry/index/detail/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/geometry/strategies/cartesian/side_by_triangle.hpp>
#include "boost_geo_serialize.hpp"

#include "hash_functions.hpp"
#include "cached_function/memoization.hpp"

void write_flattened(
        const std::string& filename,
        myriaworld::triangle_graph& g){
    using myriaworld::polar2_point;
    using myriaworld::cart3_point;
    using myriaworld::triangle_graph;
    using myriaworld::vertex_pos_t;
    using myriaworld::shared_edge_t;
    using myriaworld::country_bits_t;
    using myriaworld::vertex_fracfilled_t;
    using boost::edge_weight_t;
    using boost::edge_weight;

    boost::property_map<triangle_graph, vertex_pos_t>::type pos_map = get(vertex_pos_t(), g);
    //boost::property_map<triangle_graph, edge_weight_t>::type weight_map = get(edge_weight, g);
    //boost::property_map<triangle_graph, shared_edge_t>::type se_map = get(shared_edge_t(), g);
    boost::property_map<triangle_graph, country_bits_t>::type bs_map = get(country_bits_t(), g);
    boost::property_map<triangle_graph, vertex_fracfilled_t>::type frac_filled = get(vertex_fracfilled_t(), g);

    std::ofstream ofs(filename.c_str());
    unsigned int idx=0;
    auto vs = boost::vertices(g);
    for(auto vit = vs.first; vit != vs.second; vit++){
        {
            // do not draw triangles facing away from the viewer
            using myriaworld::cart2_point;
            namespace bg = boost::geometry;
            auto& p = pos_map[*vit].m_mappos.outer();
            cart2_point p0(p[0].get<0>(), p[0].get<1>());
            cart2_point p1(p[1].get<0>(), p[1].get<1>());
            cart2_point p2(p[2].get<0>(), p[2].get<1>());
            if(bg::strategy::side::side_by_triangle<double>::apply(p0, p1, p2) < 0)
                continue;
        }
        auto& tbits = bs_map[*vit];
        double color = frac_filled[*vit];
        for(const auto& bit : tbits){
            //color = bit.m_mapcolor;
            auto& poly = bit.m_mappos;
            if(poly.outer().size() < 2)
                continue;
            ofs << "0 " << color << " ";
            ofs << "0 " << idx++ << " ";
            for(const auto& pt : poly.outer()){
                ofs << pt.get<0>() << " " << pt.get<1>() << " ";
            }
            ofs << std::endl;
        }
    }
}

void write_s2centroids(std::string filename, myriaworld::triangle_graph& g){
    BOOST_LOG_TRIVIAL(info) << "Writing to file: " << filename;
    using myriaworld::polar2_point;
    using myriaworld::cart3_point;
    using myriaworld::triangle_graph;
    using myriaworld::vertex_pos_t;
    using myriaworld::shared_edge_t;
    using myriaworld::country_bits_t;
    using myriaworld::vertex_fracfilled_t;
    using boost::edge_weight_t;
    using boost::edge_weight;
    std::vector<polar2_point> centroids(boost::num_vertices(g));
    boost::property_map<triangle_graph, vertex_pos_t>::type pos_map = get(vertex_pos_t(), g);
    boost::property_map<triangle_graph, edge_weight_t>::type weight_map = get(edge_weight, g);
    boost::property_map<triangle_graph, shared_edge_t>::type se_map = get(shared_edge_t(), g);
    boost::property_map<triangle_graph, country_bits_t>::type bs_map = get(country_bits_t(), g);
    boost::property_map<triangle_graph, vertex_fracfilled_t>::type area_map = get(vertex_fracfilled_t(), g);
    boost::property_map<triangle_graph, vertex_fracfilled_t>::type frac_filled = get(vertex_fracfilled_t(), g);
    auto vs = boost::vertices(g);
    for(auto vit = vs.first; vit != vs.second; vit++){
        centroids[*vit] = myriaworld::geo::spherical_centroid(pos_map[*vit].m_s2_poly);
    }

    std::ofstream ofs(filename.c_str());
    auto optional_print_s2 = [&](double color, const polar2_point& a, const polar2_point& b){
        if(myriaworld::geo::crosses_dateline(a, b))
            return false;
        ofs << color << " " << a.get<0>() << " " << a.get<1>()
            << " "          << b.get<0>() << " " << b.get<1>() << std::endl;
        return true;
    };
    //auto optional_print_c2 = [&](double color, const cart3_point& a, const cart3_point& b){
        //if(a.get<2>() < 0) return;
        //if(b.get<2>() < 0) return;
        //ofs << color << " " << a.get<0>() << " " << a.get<1>()
            //<< " "          << b.get<0>() << " " << b.get<1>() << std::endl;
    //};

    if(0)
    for(auto vit = vs.first; vit != vs.second; vit++){
        //auto t = pos_map[*vit].m_mappos;
        //if(t.outer().size() != 3)
            //continue;
        //optional_print_c2(1, t.outer()[0], t.outer()[1]);
        //optional_print_c2(1, t.outer()[0], t.outer()[2]);
        //optional_print_c2(1, t.outer()[1], t.outer()[2]);
        auto t = pos_map[*vit].m_s2_poly;
        auto c = area_map[*vit];
        optional_print_s2(c, t.outer()[0], t.outer()[1]);
        optional_print_s2(c, t.outer()[0], t.outer()[2]);
        optional_print_s2(c, t.outer()[1], t.outer()[2]);
    }
    std::ofstream triapoly("triagrid_poly.txt");
    for(auto vit = vs.first; vit != vs.second; vit++){
        for(const auto& bit : bs_map[*vit]){
            const auto& b = bit.m_s2_poly.outer();
            bool ok = true;
            for (unsigned int i = 0; i < b.size(); ++i){
                ok &= optional_print_s2(10, b[i], b[(i + 1)%b.size()]);
            }
            if(!ok)
                continue;

            auto& poly = bit.m_s2_poly;
            if(poly.outer().size() < 2)
                continue;
            double color = frac_filled[*vit];
            triapoly << "0 " << color << " ";
            triapoly << "0 " << 0 << " ";
            for(const auto& pt : poly.outer()){
                triapoly << pt.get<0>() << " " << pt.get<1>() << " ";
            }
            triapoly << std::endl;
        }
    }
    auto es = boost::edges(g);
    if(1)
    for(auto eit = es.first; eit != es.second; eit++){
        auto sourcev = boost::source(*eit, g);
        //auto targetv = boost::target(*eit, g);
        //auto se = myriaworld::geo::determine_shared_vertices(pos_map[sourcev].m_s2_poly, pos_map[targetv].m_s2_poly);
        auto s = centroids[boost::source(*eit, g)];
        auto t = centroids[boost::target(*eit, g)];
        auto p = pos_map[sourcev];
        //if(se_map[*eit].is_cut)
            //optional_print_s2(1 + weight_map[*eit], s, t);
        //if(p.m_mappos.outer().size() != 3){
        //    BOOST_LOG_TRIVIAL(info) << "wrong number of vertices: " << p.m_mappos.outer().size();
        //    continue;
        //}
        //auto s = p.m_s2_poly.outer()[se.shared_l0];
        //auto t = p.m_s2_poly.outer()[se.shared_l1];
        //auto s = p.m_mappos.outer()[se.shared_l0];
        //auto t = p.m_mappos.outer()[se.shared_l1];
        if(!se_map[*eit].is_cut){
            optional_print_s2(weight_map[*eit], s, t);
        } else{
            ;//optional_print_s2(weight_map[*eit], s, t);
        }
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
myriaworld::triangle_graph side_effect_wrap2(
        std::function<void(myriaworld::triangle_graph&)> func,
        myriaworld::triangle_graph& g){
    func(g);
    return g;
}

int main(int argc, char** argv){

    namespace po = boost::program_options;
    po::options_description desc("Allowed Options");
    double sigma, wlat, wlon, clat, clon, alpha;
    std::string output, format, cmap;
    int depth, steps, rsteps;
    desc.add_options()
        ("help", "Produce help message")
        ("render", po::value<std::string>(&format)->default_value("png"), "render to this format (png, pdf, ...)")
        ("output", po::value<std::string>(&output)->default_value("myriaworld"), "an identifier for output files")
        ("cmap", po::value<std::string>(&cmap)->default_value("Paired"), "colormap (binary, ...)")
        ("unfold", po::value<int>(&steps)->default_value(1), "steps from 0 to angle")
        ("rotate", po::value<int>(&rsteps)->default_value(1), "rotation steps from (clon-180) to clon")
        ("depth", po::value<int>(&depth)->default_value(5), "How often to split triangles recursively")
        ("sigma", po::value<double>(&sigma)->default_value(0.7), "how much to smooth landmass")
        ("wlat", po::value<double>(&wlat)->default_value(0.1), "latitude graticular weight")
        ("wlon", po::value<double>(&wlon)->default_value(0.5), "longitude graticular weight")
        ("clat", po::value<double>(&clat)->default_value(10), "latitude map center")
        ("clon", po::value<double>(&clon)->default_value(10), "longitude map center")
        ("alpha", po::value<double>(&alpha)->default_value(1.), "flattening angle")
        ;
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(desc).run(), vm);
    po::notify(vm);
    
    if(vm.count("help")){
        std::cout << desc << std::endl;
        return 0;
    }

    using namespace myriaworld;
    using namespace memoization;

    disk c((boost::format("out/%s")%output).str());
    auto g = get_triangle_graph(depth);

    using std::string;
    auto countries = c(string("read_countries"), read_countries, string("easy2.txt"));
    g = c(string("country2tria"), side_effect_wrap, country2tria_s2, g, countries);
    g = c(string("smooth_landmass"), smooth_landmass, g, sigma);
    auto render = [&](int cnt){
        if(vm.count("render")){
            auto flattened_fn = boost::format("out/%s-flattened.txt") % output;
            write_flattened(flattened_fn.str(), g);
            auto output_fn = boost::format("out/%s-%04d.%s") % output % cnt % format;
            int ret = system(boost::str(boost::format("python show.py %s %s %s") % flattened_fn % output_fn % cmap).c_str());
            if(ret != 0)
                throw std::runtime_error("failed to call python renderer");
        }
    };

    int cnt=0;
    for(int step = 0; step < rsteps; step++, cnt++){
        double clon2 = clon - 180 + (step+1.0)/rsteps * 180.;
        g = c(string("determine_edge_weights"), determine_edge_weights, g, wlat, wlon, clat, clon2);
        g = c(string("determine_cuttings"), determine_cuttings, g);

        flatten(g, 0, clat, clon2);
        render(cnt);
    }

    for(int step = (steps==1?0:-1); step < steps; step++, cnt++){
        double a = (step+1.0)/steps * alpha;
        flatten(g, a, clat, clon);
        render(cnt);
    }

}
