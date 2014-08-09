#include <boost/geometry/index/rtree.hpp>
#include <boost/log/trivial.hpp>
#include <boost/progress.hpp>
#include <s2.h>
#include <s2polygon.h>
#include <s2cap.h>
#include <s2latlng.h>
#include <s2regioncoverer.h>
#include <s2polygonbuilder.h>
#include "myriaworld.h"
#include "country2tria.hpp"

namespace myriaworld
{
    void country2tria_s2(triangle_graph& g, const std::vector<country>& countries){
        BOOST_LOG_TRIVIAL(info) << "Determining overlap btw triangles, countries via Google S2";

        typedef std::shared_ptr<S2Polygon> poly_ptr;
        std::vector<poly_ptr> flat_countries;
        std::vector<std::vector<S2CellId> > flat_cells;

        BOOST_LOG_TRIVIAL(info) << "...convert from boost to S2 format";
        boost::progress_display show_progress0(countries.size());
        std::size_t n_cells = 0;
        for(const auto& c : countries){
            ++show_progress0;
            std::vector<S2Loop*> loops;
            S2PolygonBuilderOptions pbo = S2PolygonBuilderOptions::UNDIRECTED_XOR();
            S2PolygonBuilder pb(pbo);
            for(const auto& poly : c.m_s2_polys){
                std::vector<S2Point> pts;
                for(const auto& v : poly.outer()){
                    S2LatLng s2p(S2LatLng::FromDegrees(
                                v.get<0>(),
                                v.get<1>()));
                    //assert(s2p.is_valid());
                    s2p = s2p.Normalized();
                    pts.push_back(s2p.ToPoint());
                }
                //pts.push_back(pts[0]);
                //assert(S2::RobustCCW(pts[0], pts[1], pts[2]) != 0);
                S2Loop loop(pts);
                pb.AddLoop(&loop);
                //loops.push_back(new S2Loop(pts));
            }
            //auto cptr = std::make_shared<S2Polygon>(&loops);
            auto cptr = std::make_shared<S2Polygon>();
            pb.AssemblePolygon(cptr.get(), NULL);
            flat_countries.push_back(cptr);

            S2RegionCoverer coverer;
            double diameter = 2*cptr->GetCapBound().angle().radians();
            int min_level = S2::kMaxWidth.GetMinLevel(diameter);
            coverer.set_min_level(min_level);
            coverer.set_max_level(min_level+4);
            coverer.set_max_cells(50);

            std::vector<S2CellId> cells;
            coverer.GetCovering(*cptr, &cells);
            flat_cells.push_back(cells);
            n_cells += cells.size();
        }
        BOOST_LOG_TRIVIAL(info) << "avg cells/country: "<<n_cells/(float)countries.size();
        boost::property_map<triangle_graph, vertex_pos_t>::type pos_map = get(vertex_pos_t(), g);
        boost::property_map<triangle_graph, country_bits_t>::type bs_map = get(country_bits_t(), g);
        BOOST_LOG_TRIVIAL(info) << "...intersection operation...";
        auto vs = boost::vertices(g);
        boost::progress_display show_progress(boost::num_vertices(g));
        for (auto vit = vs.first; vit != vs.second; ++vit)
        {
            ++show_progress;

            S2Polygon tria;
            {
                std::vector<S2Point> v;
                for (const auto& pt : pos_map[*vit].m_s2_poly.outer())
                {
                    S2LatLng s2p(S2LatLng::FromDegrees(
                                pt.get<1>(),
                                pt.get<0>()));
                    s2p = s2p.Normalized();
                    v.push_back(s2p.ToPoint());
                }
                if(S2::SimpleCCW(v[0], v[1], v[2])<0)
                    std::swap(v[0], v[2]);
                S2Loop* loop = new S2Loop(v);
                std::vector<S2Loop*> loops;
                loops.push_back(loop);
                tria.Init(&loops);
            }

            std::vector<S2CellId> triacells;
            {
                S2RegionCoverer coverer;
                double diameter = tria.GetCapBound().angle().radians();
                int min_level = S2::kMaxWidth.GetMinLevel(diameter);
                coverer.set_min_level(min_level);
                coverer.set_max_level(min_level+5);
                coverer.set_max_cells(5);
                coverer.GetCovering(tria, &triacells);
            }
            //for(const auto& cptr : flat_countries){
            for (unsigned int cidx = 0; cidx < flat_countries.size(); ++cidx)
            {
                std::shared_ptr<S2Polygon> cptr = flat_countries[cidx];

                bool mayintersect = false;
                for(const auto& id0 : flat_cells[cidx]){
                    S2Cell cell0(id0);
                    S2Polygon cellp0(cell0);
                    for(const auto& id1 : triacells){
                        S2Cell cell1(id1);
                        if(cellp0.MayIntersect(cell1)){
                            mayintersect = true;
                            break;
                        }
                    }
                    if(mayintersect)
                        break;
                }
                if(!mayintersect)
                    continue;

#if 0
                S2Polygon isect;
                isect.InitToIntersectionSloppy(cptr.get(), &tria, S1Angle::Degrees(0.000001));
                //isect.InitToIntersection(cptr.get(), &tria);
#else
                std::vector<S2Polygon*> pieces;
                for (const auto& id : flat_cells[cidx])
                {
                    S2Cell cell(id);
                    if(tria.MayIntersect(cell)){
                        S2Polygon cpiece;
                        S2Polygon cellpoly(cell);
                        cpiece.InitToIntersection(cptr.get(), &cellpoly);
                        if(cpiece.num_loops() == 0)
                            continue;
                        auto piece = new S2Polygon();
                        piece->InitToIntersection(&cpiece, &tria);
                        if(piece->num_loops() == 0)
                            delete piece;
                        else
                            pieces.push_back(piece);
                    }
                }
                std::shared_ptr<S2Polygon> isect_ptr( S2Polygon::DestructiveUnion(&pieces) );
                S2Polygon& isect = *isect_ptr;
#endif
                int n_loops = isect.num_loops();
                if(n_loops == 0)
                    continue;
                for (int i = 0; i < n_loops; ++i)
                {
                    myriaworld::polar2_polygon s2poly;
                    S2Loop& loop = *isect.loop(i);
                    //if(loop.is_hole())
                        //continue;
                    //assert(!loop.is_hole());
                    int n_vert = loop.num_vertices();
                    // our polygons are not closed!
                    for (int j = 0; j < n_vert; ++j)
                    {
                        S2Point p = loop.vertex(j);
                        S2LatLng s2p(p);
                        s2poly.outer().push_back(
                                polar2_point(s2p.lng().degrees(),
                                    s2p.lat().degrees()));
                    }
                    country_bit bit;
                    bit.m_s2_poly = s2poly;
                    bs_map[*vit].push_back(bit);
                }
            }
        }
        BOOST_LOG_TRIVIAL(info) << "...done.";
    }
}
