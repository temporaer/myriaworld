#include "myriaworld.h"
#include <iostream>
#include <fstream>
#include <limits>
#include <boost/progress.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/geometry/strategies/transform.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/log/trivial.hpp>
#include <boost/algorithm/string/split.hpp>
#include "geo.hpp"
#include "SPHEmesh.hpp"


namespace bg = boost::geometry;

namespace myriaworld{
    namespace detail{
        struct max_depth_reached{};
        class gauss_accu : public boost::default_bfs_visitor
        {
            public:
                gauss_accu(polar2_point p, double* sum, double* wsum, double* maxd, int* n_nodes, double sigma):m_cpos(p),m_weightsum(wsum),m_sum(sum),m_sigma(sigma),m_maxd(maxd),m_n_nodes(n_nodes){}
                polar2_point m_cpos;
                double *m_weightsum, *m_sum;
                double m_sigma;
                double* m_maxd;
                int* m_n_nodes;
                template < typename Vertex, typename Graph >
                    void examine_vertex(Vertex u, const Graph & g)
                    {
                        //boost::property_map<triangle_graph, vertex_pos_t>::const_type pos_map 
                        //    = get(vertex_pos_t(), g);
                        boost::property_map<triangle_graph, vertex_centroid_t>::const_type centroid_map 
                            = get(vertex_centroid_t(), g);
                        boost::property_map<triangle_graph, vertex_fracfilled_t>::const_type
                            frac_filled = get(vertex_fracfilled_t(), g);
                        boost::property_map<triangle_graph, vertex_area_t>::const_type area_map 
                            = get(vertex_area_t(), g);
                        const auto& p0 = centroid_map[u];
                        const auto& p1 = m_cpos;
                        double dist = geo::haversine_distance(p0, p1);
                        if(dist > M_PI/2.)
                            throw max_depth_reached();
                        (*m_n_nodes) ++;
                        double weight = exp(-dist*dist / (m_sigma * m_sigma));
                        if(weight < 1e-3)
                            throw max_depth_reached();
                        *m_maxd = std::min(weight, *m_maxd);
                        //BOOST_LOG_TRIVIAL(debug) << "        w: " << weight;
                        *m_sum += weight * area_map[u] * frac_filled[u];
                        *m_weightsum += weight * area_map[u];
                    }
        };
    }
}


namespace myriaworld
{
    std::vector<country> read_countries(const std::string& filename){
        std::ifstream ifs(filename.c_str());
        std::vector<country> ret;
        //unsigned int cnt=0;
        while(ifs){
            std::string line;
            std::vector<std::string> v;
            std::getline(ifs, line);
            boost::split(v, line, boost::is_any_of(";"));
            if(v.size() != 4)
                break;

            country c;
            c.m_name = v[0];
            c.m_mapcolor = boost::lexical_cast<double>(v[1]);
            c.m_scalerank = boost::lexical_cast<double>(v[2]);

            BOOST_LOG_TRIVIAL(debug) << "Reading " << c.m_name;
            if(v[3].find("MULTIPOLY") != std::string::npos){
                bg::read_wkt(v[3], c.m_s2_polys);
                for(auto & p : c.m_s2_polys){
                    // wkt is closed, we want non-closed entries
                    p.outer().pop_back();
                    // WKT is CW, Google S2 wants CCW.
                    std::reverse(p.outer().begin(), p.outer().end());
                }
            }else{
                polar2_polygon poly2;
                bg::read_wkt(v[3], poly2);

                // wkt is closed, we want non-closed entries
                poly2.outer().pop_back();
                // WKT is CW, Google S2 wants CCW.
                std::reverse(poly2.outer().begin(), poly2.outer().end());

                c.m_s2_polys.push_back(poly2);
            }
            ret.push_back(c);
            //if(cnt++ > 10) // TODO remove
                //break;
        }
        return ret;
    }
    std::vector<country> read_countries2(const std::string& filename){
        std::ifstream ifs(filename.c_str());
        unsigned n_poly, n_pts;
        double lat, lng;
        double minlat= std::numeric_limits<double>::infinity(),
               maxlat=-std::numeric_limits<double>::infinity();
        double minlng= std::numeric_limits<double>::infinity(),
               maxlng=-std::numeric_limits<double>::infinity();
        std::vector<country> countries;
        int cidx = 0;
        while(ifs){
            ifs >> n_poly;
            if(n_poly == 0)
                continue;
            country c;
            for(unsigned r=0; r< n_poly; r++){
                ifs >> n_pts;
                polar2_polygon poly;
                for(unsigned p=0; p<n_pts; p++){
                    ifs >> lng >> lat;
                    minlat = std::min(minlat, lat);
                    minlng = std::min(minlng, lng);
                    maxlat = std::max(maxlat, lat);
                    maxlng = std::max(maxlng, lng);
                    append(poly, polar2_point(lat, lng));
                }

                // NOTE we assume that all countries are specified
                // in ccw order

                //polar2_polygon poly2;
                //boost::geometry::simplify(poly, poly2, 0.001);
                //boost::geometry::correct(poly2);
                //if(poly2.outer().size() > 4)
                    //c.m_s2_polys.push_back(poly2);
                c.m_s2_polys.push_back(poly);
            }
            //assert(c.size() == n_poly);
            //BOOST_LOG_TRIVIAL(debug) << "# country...";
            //if(cidx <= 36)
            if(c.m_s2_polys.size() > 0)
                countries.push_back(c);
            //assert(countries.back().size() == n_poly);
            cidx ++;
        }
        BOOST_LOG_TRIVIAL(info) << "countries -- minlat: "<< minlat << " maxlat: "<<maxlat;
        BOOST_LOG_TRIVIAL(info) << "countries -- minlng: "<< minlng << " maxlng: "<<maxlng;
        return countries;
    }

    std::vector<polar2_polygon> read_triangle_grid(const std::string& filename){
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

            {
                boost::geometry::strategy::transform::from_cartesian_3_to_spherical_equatorial_2<cart3_point, polar2_point> strategy;
                for (unsigned int i = 0; i < triangles.size(); ++i)
                {
                    boost::geometry::transform(ctriangles[i], triangles[i], strategy);

                    //BOOST_LOG_TRIVIAL(info) << "Triangle Area: " << bg::area(ctriangles[i]);
                    //BOOST_LOG_TRIVIAL(info) << "Triangle Area:" << bg::area(triangles[i]);
                    //triangles[i] = geo::rotated(triangles[i], 2.1828, 2.1725);
                }
            }

            triangle_graph tg(triangles.size());

            //boost::property_map<triangle_graph, vertex_index_t>::type index_map 
            //= get(vertex_index, tg);

            boost::property_map<triangle_graph, vertex_pos_t>::type pos_map 
                = get(vertex_pos_t(), tg);

            boost::property_map<triangle_graph, n_shared_vert_t>::type n_shared_vert_map 
                = get(n_shared_vert_t(), tg);

            // the position, the vertex ID, the triangle ID
            typedef boost::tuple<cart2_point, unsigned, unsigned> value;
            bgi::rtree <value, bgi::quadratic<16> > rtree;
            std::vector<polar2_point> vertices;
            std::vector<std::vector<unsigned int> > 
                vertex_to_tria(boost::num_vertices(tg));

            auto add_node_checked = [&](unsigned int tria, const polar2_point& a){
                std::vector<value> res;
                using namespace myriaworld;

                cart2_point ca;
                geo::cast_polar2_equatorial_2_cart_2(ca, a, 0.);

                rtree.query(bgi::nearest(ca, 1), std::back_inserter(res));
                if(res.size() == 0 || geo::haversine_distance(vertices[res[0].get<1>()], a) > 0.001){
                    unsigned int idx = vertices.size();
                    vertices.push_back(a);
                    rtree.insert({ca, idx, tria});
                    vertex_to_tria[idx] = {tria};
                }else{
                    // we found another use of this vertex.
                    unsigned int other = res.front().get<1>();
                    for(auto tria2 : vertex_to_tria[other]){
                        auto e = boost::edge(tria2, tria, tg);
                        if(!e.second){
                            boost::add_edge(tria2, tria, tg);
                            n_shared_vert_map[boost::edge(tria2, tria, tg).first] = 1;
                        }else{
                            n_shared_vert_map[e.first]++;
                            assert( n_shared_vert_map[e.first] <= 2);
                        }
                    }
                    vertex_to_tria[other].push_back(tria);
                }

            };

            unsigned int tria_idx = 0;
            for(const auto& c : triangles){
                pos_map[tria_idx].m_s2_poly.outer().push_back(c.outer()[0]);
                pos_map[tria_idx].m_s2_poly.outer().push_back(c.outer()[1]);
                pos_map[tria_idx].m_s2_poly.outer().push_back(c.outer()[2]);

                assert(c.outer().size() == 3);
                add_node_checked(tria_idx, c.outer()[0]);
                add_node_checked(tria_idx, c.outer()[1]);
                add_node_checked(tria_idx, c.outer()[2]);
                ++tria_idx;
            }

            // remove all edges which have n_shared_vert != 2
            boost::remove_edge_if(
                    [&](const triaedge_descriptor& e){
                        return n_shared_vert_map[e] < 2;
                    }, tg);
            return tg;
        }

    myriaworld::triangle_graph
        smooth_landmass(myriaworld::triangle_graph& g, double sigma){
            using namespace myriaworld;
            namespace bg = boost::geometry;
            using boost::vertex_index_t;
            using boost::vertex_index;
            using boost::edge_weight_t;
            using boost::edge_weight;
            using boost::property_map;
            using boost::vertices;
            using boost::edges;


            //property_map<triangle_graph, vertex_index_t>::type vertex_index_map = get(vertex_index, g);
            //property_map<triangle_graph, vertex_pos_t>::type pos_map = get(vertex_pos_t(), g);
            property_map<triangle_graph, vertex_centroid_t>::type centroid_map = get(vertex_centroid_t(), g);
            //property_map<triangle_graph, vertex_area_t>::type area_map = get(vertex_area_t(), g);
            property_map<triangle_graph, vertex_fracfilled_t>::type frac_filled = get(vertex_fracfilled_t(), g);
            property_map<triangle_graph, edge_weight_t>::type weightmap = get(edge_weight, g);
            //property_map<triangle_graph, country_bits_t>::type bits_map = get(country_bits_t(), g);

            // initialize weights to constant for BFS
            auto es = edges(g);
            for(auto eit = es.first; eit != es.second; eit++){
                weightmap[*eit] = 1.0;
            }

            // 2. smooth triangle filled-counters
            if(true){
                BOOST_LOG_TRIVIAL(info) << "Smooth Landmass";
                auto vs = vertices(g);
                boost::progress_display show_progress(std::distance(vs.first, vs.second));
                std::vector<double> frac_filled2(boost::num_vertices(g));
                for(auto vit = vs.first; vit!=vs.second; vit++){
                    ++show_progress;
                    //BOOST_LOG_TRIVIAL(debug) << "FF before: " << frac_filled[*vit];
                    double sum=0, wsum=0, maxd=1E6;
                    int n_nodes=0;
                    detail::gauss_accu acc(centroid_map[*vit], &sum, &wsum, &maxd, &n_nodes, sigma);
                    try{
                        breadth_first_search(g, *vit, visitor(acc));
                    }catch(detail::max_depth_reached){
                        ;
                    }
                    frac_filled2[*vit] = std::max(0., sum / (wsum + 0.000001));
                    static int cnt=0;
                    if(cnt++ % 100 == 0){
                        //BOOST_LOG_TRIVIAL(debug) << "  maxd : " << maxd << " n_nodes: " << n_nodes;
                    }
                }
                for(unsigned int i=0; i < boost::num_vertices(g); i++)
                    frac_filled[i] = frac_filled2[i];
            }
            return g;
        }

    myriaworld::triangle_graph
        determine_edge_weights(myriaworld::triangle_graph& g, double wlat, double wlon, double clat, double clon){
            using namespace boost;
            property_map<triangle_graph, vertex_index_t>::type vertex_index_map = get(vertex_index, g);
            property_map<triangle_graph, vertex_fracfilled_t>::type frac_filled = get(vertex_fracfilled_t(), g);
            property_map<triangle_graph, vertex_centroid_t>::type centroid_map = get(vertex_centroid_t(), g);
            property_map<triangle_graph, edge_weight_t>::type weightmap = get(edge_weight, g);
            property_map<triangle_graph, vertex_area_t>::type area_map = get(vertex_area_t(), g);

            if(0){
                // strategically increase the weight of some countries,
                // which are notorious for being split
                property_map<triangle_graph, country_bits_t>::type bits_map = get(country_bits_t(), g);
                // main work is done in country2tria_s2.cpp
                auto vs = vertices(g);
                for(auto vit = vs.first; vit!=vs.second; vit++){
                    bool double_weight = false;
                    for(const auto& p_ : bits_map[*vit]){
                        if(p_.m_name == "Greenland"){
                            double_weight = true;
                            break;
                        }
                        if(p_.m_name == "Antarctica"){
                            double_weight = true;
                            break;
                        }
                    }
                    if(double_weight)
                        frac_filled[*vit] *= 2.;
                }
            }

            auto W0 = [&](double v, double v0){
                return pow(fabs(v - v0)/180., 0.5);
            };
            auto W1 = [&](double v, double v0){
                return pow(fabs(v - v0)/45., 0.3);
            };
            auto WNorm = [&](double a, double b){
                return wlat * a*a + wlon * b*b;
            };
            {
                // 3. set the edge weights.
                BOOST_LOG_TRIVIAL(info) << "Set edge weights";
                auto es = edges(g);
                for(auto eit = es.first; eit != es.second; eit++){
                    const auto sourcev = source(*eit, g);
                    const auto targetv = target(*eit, g);
                    double sum = 
                        frac_filled[vertex_index_map[sourcev]] * area_map[sourcev] +
                        frac_filled[vertex_index_map[targetv]] * area_map[targetv];
                    double weightsum = area_map[sourcev] + area_map[targetv];
                    sum /= weightsum; 

                    sum = (1 - sum) * (1 - sum)
                        * 
                        WNorm( W0(centroid_map[sourcev].get<0>(), clat)
                          , W1(centroid_map[sourcev].get<1>(), clon));
                    //weightmap[*eit] = std::max(0., sum);
                    weightmap[*eit] = std::exp(sum);
                }
            }

            return g;
    }
}
