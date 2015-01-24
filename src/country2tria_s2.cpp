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
//#include <szl/hashutils.hpp>
//#include <strutil.hpp>

namespace myriaworld
{
    void country2tria_s2(triangle_graph& g, const std::vector<country>& countries){
        BOOST_LOG_TRIVIAL(info) << "Determining overlap btw triangles, countries via Google S2";

        typedef std::shared_ptr<S2Polygon> poly_ptr;
        std::vector<std::vector<poly_ptr> > flat_countries;
        std::vector<std::vector<S2CellId> > flat_cells;

        BOOST_LOG_TRIVIAL(info) << "...convert from boost to S2 format";
        boost::progress_display show_progress0(countries.size());
        std::size_t n_cells = 0;
        for(const auto& c : countries){
            //++show_progress0;
            std::vector<S2Polygon*> polys;
            BOOST_LOG_TRIVIAL(info) << c.m_name;
            for(const auto& poly : c.m_s2_polys){
                std::vector<S2Point> pts;
                for(const auto& v : poly.outer()){
                    S2LatLng s2p(S2LatLng::FromDegrees(
                                v.get<1>(),
                                v.get<0>()));
                    s2p = s2p.Normalized();
                    pts.push_back(s2p.ToPoint());
                }

                S2Loop* loop = new S2Loop(pts);
                if(loop->GetArea() > 10){
                    delete loop;
                    std::reverse(pts.begin(), pts.end());
                    loop = new S2Loop(pts);
                }
                std::vector<S2Loop*> lpv(1, loop);
                polys.push_back(new S2Polygon(&lpv));
            }
            S2Polygon* tmp = S2Polygon::DestructiveUnion(&polys);
            auto cptr = std::shared_ptr<S2Polygon>(tmp);
            double carea = cptr->GetArea();
            //BOOST_LOG_TRIVIAL(info) << "Country Area: " << carea;

            // now take the country apart into cells covering it
            S2RegionCoverer coverer;
            double diameter = 2*cptr->GetCapBound().angle().radians();
            int min_level = S2::kMaxWidth.GetMinLevel(diameter);
            coverer.set_min_level(min_level);
            coverer.set_max_level(min_level+3);
            coverer.set_max_cells(50);

            std::vector<S2CellId> cells;
            coverer.GetCovering(*cptr, &cells);

            // For each cell intersecting the country, determine the piece
            // which is in that cell.
            std::vector<poly_ptr> flat_country;
            std::vector<S2CellId> flat_cell;
            for(const auto& id0 : cells){
                auto piece = std::make_shared<S2Polygon>();
                S2Cell cell0(id0);
                S2Polygon cellp0(cell0);
                piece->InitToIntersection(cptr.get(), &cellp0);
                assert(piece->GetArea() < 10);
                //BOOST_LOG_TRIVIAL(info) << "Country Area: " << carea;
                //BOOST_LOG_TRIVIAL(info) << "Piece Area: " << piece->GetArea();
                // the following two resulted in small holes around countries
                //if(piece->GetArea() < 0.00001)
                //    continue;
                //if(piece->GetArea() > carea + 0.0000001)
                //    continue;
                flat_country.push_back(piece);
                flat_cell.push_back(id0);
            }
            flat_countries.push_back(flat_country);
            flat_cells.push_back(flat_cell);
            n_cells += cells.size();
        }
        //BOOST_LOG_TRIVIAL(info) << "avg cells/country: "<<n_cells/(float)countries.size();
        boost::property_map<triangle_graph, vertex_pos_t>::type pos_map = get(vertex_pos_t(), g);
        boost::property_map<triangle_graph, country_bits_t>::type bs_map = get(country_bits_t(), g);
        boost::property_map<triangle_graph, vertex_centroid_t>::type centroid_map = get(vertex_centroid_t(), g);
        boost::property_map<triangle_graph, vertex_fracfilled_t>::type frac_filled = get(vertex_fracfilled_t(), g);
        boost::property_map<triangle_graph, vertex_area_t>::type area_map = get(vertex_area_t(), g);
        //BOOST_LOG_TRIVIAL(info) << "...intersection operation...";
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
                //if(S2::SimpleCCW(v[0], v[1], v[2])<0)
                    //std::swap(v[0], v[2]);
                S2Loop* loop = new S2Loop(v);
                assert(loop->GetArea() < 0.1);
                std::vector<S2Loop*> loops;
                loops.push_back(loop);
                tria.Init(&loops);
            }
            auto centroid = tria.GetCentroid();
            centroid_map[*vit] = polar2_point(
                                S2LatLng::Longitude(centroid).degrees(),
                                S2LatLng::Latitude(centroid).degrees());
            double tria_area = tria.GetArea();
            double tria_country_area = 0.;

            std::vector<S2Polygon*> allbits;
            for (unsigned int cidx = 0; cidx < flat_countries.size(); ++cidx)
            {
#if 0
                S2Polygon isect;
                //isect.InitToIntersectionSloppy(cptr.get(), &tria, S1Angle::Degrees(0.000001));
                isect.InitToIntersection(cptr.get(), &tria);
#else
                std::vector<S2Polygon*> pieces;
                for (unsigned int cellidx = 0; cellidx < flat_countries[cidx].size(); ++cellidx)
                {
                    std::shared_ptr<S2Polygon> cptr = flat_countries[cidx][cellidx];
                    const auto& id = flat_cells[cidx][cellidx];
                    
                    S2Cell cell(id);
                    if(tria.MayIntersect(cell)){
                        auto piece = new S2Polygon();
                        piece->InitToIntersection(cptr.get(), &tria);
                        //BOOST_LOG_TRIVIAL(info) << "---------------------";
                        //BOOST_LOG_TRIVIAL(info) << "Country     : " << cptr->GetArea();
                        //BOOST_LOG_TRIVIAL(info) << "Tria        : " << tria.GetArea();
                        //BOOST_LOG_TRIVIAL(info) << "Intersection: " << piece->GetArea();
                        assert(piece->GetArea() < 10);
                        assert(piece->GetArea() < cptr->GetArea() + tria.GetArea());
                        if(piece->num_loops() == 0)
                            delete piece;
                        else
                            pieces.push_back(piece);
                    }
                }
                S2Polygon* isect = S2Polygon::DestructiveUnion(&pieces);
#endif
                int n_loops = isect->num_loops();
                if(n_loops == 0)
                    continue;
                //if(isect->GetArea() < 0.00001)
                //    continue;
                //BOOST_LOG_TRIVIAL(info) << "Union:        " << isect->GetArea();
                assert(isect->GetArea() < 10);
                assert(isect->GetArea() <= tria.GetArea() + 0.000001);
                tria_country_area += isect->GetArea();
                allbits.push_back(isect);
                for (int i = 0; i < n_loops; ++i)
                {
                    myriaworld::polar2_polygon s2poly;
                    S2Loop& loop = *isect->loop(i);
                    if(loop.is_hole())
                        continue;
                    //assert(!loop.is_hole());
                    int n_vert = loop.num_vertices();
                    // our polygons are not closed!
                    for (int j = 0; j < n_vert; ++j)
                    {
                        S2Point p = loop.vertex(j);
                        S2LatLng s2p(p);
                        s2poly.outer().push_back(
                                polar2_point(
                                    s2p.lng().degrees(),
                                    s2p.lat().degrees()));
                    }
                    if(n_vert < 3) continue;
                    country_bit bit(countries[cidx]); // copies color infos
                    //bg::correct(s2poly);
                    bit.m_s2_poly = s2poly;
                    bs_map[*vit].push_back(bit);
                }
            }

            area_map[*vit] = tria_area;
            double frac = tria_country_area / (tria_area + 0.0000001);
            if(frac < 0 || frac > 1.001)
                BOOST_LOG_TRIVIAL(info) << "frac filled: " << frac;
            frac_filled[*vit] = std::max(0.0, std::min(1.0, tria_country_area / (tria_area + 0.0000001)));

            auto allbits_ = std::shared_ptr<S2Polygon>(
                    S2Polygon::DestructiveUnion(&allbits));

            // everything else must be water!
            S2Polygon waterbits;
            waterbits.InitToDifference(&tria, allbits_.get());
            for (int i = 0; i < waterbits.num_loops(); ++i)
            {
                myriaworld::polar2_polygon s2poly;
                S2Loop& loop = *waterbits.loop(i);
                if(loop.is_hole())
                    continue;
                //assert(!loop.is_hole());
                int n_vert = loop.num_vertices();
                // our polygons are not closed!
                for (int j = 0; j < n_vert; ++j)
                {
                    S2Point p = loop.vertex(j);
                    S2LatLng s2p(p);
                    s2poly.outer().push_back(
                            polar2_point(
                                s2p.lng().degrees(),
                                s2p.lat().degrees()));
                }
                if(n_vert < 3) continue;
                country_bit bit;
                bit.m_scalerank = 0;
                bit.m_mapcolor = -1;
                bit.m_name = "OCEAN";
                //bg::correct(s2poly);
                bit.m_s2_poly = s2poly;
                bs_map[*vit].push_back(bit);
            }
        }
        BOOST_LOG_TRIVIAL(info) << "...done.";
    }
}
