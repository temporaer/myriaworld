#include "myriaworld.h"
#include "geo.hpp"
#include <boost/numeric/bindings/traits/ublas_matrix.hpp> 
#include <boost/math/quaternion.hpp>
#include <boost/log/trivial.hpp>

// this is from boost numeric bindings!
#include <boost/numeric/bindings/atlas/clapack.hpp>
#include<boost/geometry/strategies/transform/matrix_transformers.hpp>


namespace myriaworld
{
    namespace geo
    {
        bool crosses_dateline(const polar2_point& a, const polar2_point& b){
            bool crosses = false;
            double lat0 = a.get<0>();
            double lat1 = b.get<0>();
            if(lat0 * lat1 < 0
                    && fabs(lat0) > 170 && fabs(lat1) > 170){
                crosses = true;
            }
            return crosses;
        }
        bool crosses_dateline(const polar2_polygon& poly){
            bool crosses = false;
            for(unsigned int i=0; i < poly.outer().size()-1; ++i){
                double lat0 = poly.outer()[i+0].get<0>();
                double lat1 = poly.outer()[i+1].get<0>();
                if(lat0 * lat1 < 0
                        && fabs(lat0) > 170 && fabs(lat1) > 170){
                    crosses = true;
                    break;
                }
            }
            return crosses;
        }

    void cast_polar2_equatorial_2_cart_2(cart2_point& dst, const polar2_point& src, double addlat, double addlng){
        double f = src.get<0>() + addlat;
        while(f >  180.) f -= 360.;
        while(f < -180.) f += 360.;
        dst.set<0>(f);

        f = src.get<1>() + addlng;
        while(f >  180.) f -= 360.;
        while(f < -180.) f += 360.;
        dst.set<1>(f);
    }
    void cast_polar2_equatorial_2_cart_2(cart2_polygon& dst, const polar2_polygon& src, double addlat, double addlng){
        dst.outer().clear();
        for(const auto& pt : src.outer()){
            cart2_point c2p;
            cast_polar2_equatorial_2_cart_2(c2p, pt, addlat, addlng);
            bg::append(dst.outer(), c2p);
        }
        assert(dst.outer().size() == src.outer().size());
    }
    void cast_polar2_equatorial_2_cart_2(cart2_multipolygon& dst, const polar2_multipolygon& src, double addlat, double addlng){
        dst.clear();
        for(const auto& p : src){
            cart2_polygon ep;
            cast_polar2_equatorial_2_cart_2(ep, p, addlat, addlng);
            dst.push_back(ep);
        }
        assert(dst.size() == src.size());
    }

    void cast_cart_2_polar2_equatorial_2(polar2_point& dst, const cart2_point& src, double addlat, double addlng){
        double f = src.get<0>() + addlat;
        while(f >  180.) f -= 360.;
        while(f < -180.) f += 360.;
        dst.set<0>(f);

        f = src.get<1>() + addlng;
        while(f >  180.) f -= 360.;
        while(f < -180.) f += 360.;
        dst.set<1>(f);
    }
    void cast_cart_2_polar2_equatorial_2(polar2_polygon& dst, const cart2_polygon& src, double addlat, double addlng){
        dst.outer().clear();
        for(const auto& pt : src.outer()){
            polar2_point s2p;
            cast_cart_2_polar2_equatorial_2(s2p, pt, addlat, addlng);
            bg::append(dst, s2p);
        }
    }
    void cast_cart_2_polar2_equatorial_2(polar2_multipolygon& dst, const cart2_multipolygon& src, double addlat, double addlng){
        dst.clear();
        for(const auto& p : src){
            polar2_polygon sp;
            cast_cart_2_polar2_equatorial_2(sp, p, addlat, addlng);
            dst.push_back(sp);
        }
    }
        //polar2_point rotated(const polar2_point& sp, double lat_off, double lng_off){
        ////polar2_point rotated(const polar2_point& , double , double ){
        //    using namespace boost::geometry::strategy::transform;
        //    cart3_point c3p, tmp;
        //    polar2_point sp2;
        //    from_spherical_equatorial_2_to_cartesian_3<polar2_point, cart3_point> sph2cart3;
        //    from_cartesian_3_to_spherical_equatorial_2<cart3_point, polar2_point> cart3sph2;

        //    auto quaternion = boost::math::spherical(1.0, 0., lat_off/180.*M_PI, lng_off/180.*M_PI);
        //    double qw = quaternion.R_component_1();
        //    double qx = quaternion.R_component_2();
        //    double qy = quaternion.R_component_3();
        //    double qz = quaternion.R_component_4();
        //    const double n = 1.0f/sqrt(qx*qx+qy*qy+qz*qz+qw*qw);
        //    qx *= n; qy *= n; qz *= n; qw *= n;
        //    boost::geometry::strategy::transform::ublas_transformer<double, 3, 3> combined(
        //            1.0f - 2.0f*qy*qy - 2.0f*qz*qz, 2.0f*qx*qy - 2.0f*qz*qw, 2.0f*qx*qz + 2.0f*qy*qw, 0.0f,
        //            2.0f*qx*qy + 2.0f*qz*qw, 1.0f - 2.0f*qx*qx - 2.0f*qz*qz, 2.0f*qy*qz - 2.0f*qx*qw, 0.0f,
        //            2.0f*qx*qz - 2.0f*qy*qw, 2.0f*qy*qz + 2.0f*qx*qw, 1.0f - 2.0f*qx*qx - 2.0f*qy*qy, 0.0f,
        //            0.0f, 0.0f, 0.0f, 1.0f);
            
        //    bg::transform(sp, c3p, sph2cart3);
        //    bg::transform(c3p, tmp, combined);
        //    bg::transform(tmp, sp2, cart3sph2);
        //    return sp2;
        //}
        //polar2_polygon rotated(const polar2_polygon&, double , double ){
        cart3_polygon rotated(const cart3_polygon& cp, double lat_off, double lng_off, double roll){
            using namespace boost::geometry::strategy::transform;
            cart3_polygon res;
            
            double a1 = lat_off/180.*M_PI;
            double a2 = roll/180.*M_PI;
            double a3 = lng_off/180.*M_PI;

            double c1 = cos(a1), s1 = sin(a1);
            double c2 = cos(a2), s2 = sin(a2);
            double c3 = cos(a3), s3 = sin(a3);

            boost::geometry::strategy::transform::ublas_transformer<double, 3, 3> combined(
                    c2*c3,            -c2*s3,           s2,     0.,
                    c1*s3 + c3*s1*s2, c1*c3 - s1*s2*s3, -c2*s1, 0.,
                    s1*s3 - c1*c3*s2, c3*s1 + c1*s2*s3, c1*c2,  0.,
                    0.,               0.,               0.,     1.
                    );
            
            cart3_point tmp;
            for(const auto& p : cp.outer()){
                bg::transform(p, tmp, combined);
                bg::append(res, tmp);
            }
            return res;
        }
        polar2_polygon rotated(const polar2_polygon& sp, double lat_off, double lng_off, double roll){
            using namespace boost::geometry::strategy::transform;
            from_spherical_equatorial_2_to_cartesian_3<polar2_point, cart3_point> sph2cart3;
            from_cartesian_3_to_spherical_equatorial_2<cart3_point, polar2_point> cart3sph2;

            polar2_polygon res;
            cart3_polygon c3p, tmp;
            bg::transform(sp, c3p, sph2cart3);
            tmp = rotated(c3p, lat_off, lng_off, roll);
            bg::transform(tmp, res, cart3sph2);
            return res;
        }

        double dot(const cart3_point& a, const cart3_point& b){
            return (
                    a.get<0>() * b.get<0>()+
                    a.get<1>() * b.get<1>()+
                    a.get<2>() * b.get<2>());
        }
        cart3_point difference(const polar2_point& a, const polar2_point& b){
            return cart3_point(
                    a.get<0>() - b.get<0>(),
                    a.get<1>() - b.get<1>(), 0.);
        }
        cart3_point difference(const cart3_point& a, const cart3_point& b){
            return cart3_point(
                    a.get<0>() - b.get<0>(),
                    a.get<1>() - b.get<1>(),
                    a.get<2>() - b.get<2>());
        }
        cart3_point sum(const cart3_point& a, const cart3_point& b){
            return cart3_point(
                    a.get<0>() + b.get<0>(),
                    a.get<1>() + b.get<1>(),
                    a.get<2>() + b.get<2>());
        }
        cart3_point cross(const cart3_point& a, const cart3_point& b){
            return cart3_point(
                    a.get<1>() * b.get<2>() - a.get<2>() * b.get<1>(),
                    a.get<2>() * b.get<0>() - a.get<0>() * b.get<2>(),
                    a.get<0>() * b.get<1>() - a.get<1>() * b.get<0>()
                    );
        }
        cart3_point negate(const cart3_point& a){
            return cart3_point( -a.get<0>(), -a.get<1>(), -a.get<2>());
        }
        cart3_point divide(const cart3_point& a, double div){
            return cart3_point( 
                    a.get<0>()/div,
                    a.get<1>()/div,
                    a.get<2>()/div);
        }
        double norm2(const cart3_point& a){
            return (  a.get<0>() * a.get<0>() 
                    + a.get<1>() * a.get<1>() 
                    + a.get<2>() * a.get<2>());
        }
        double norm(const cart3_point& a){
            return std::sqrt(norm2(a));
        }
        double  angle_between(
                const cart3_polygon& t0, 
                const cart3_polygon& t1){
            using namespace boost::geometry::strategy::transform;
            cart3_point n0 = cross(
                    difference(t0.outer()[1], t0.outer()[0]), 
                    difference(t0.outer()[2], t0.outer()[0]));
            cart3_point n1 = cross(
                    difference(t1.outer()[1], t1.outer()[0]), 
                    difference(t1.outer()[2], t1.outer()[0]));

            //return atan2(norm(cross(n0, n1)), dot(n0, n1));

            //if(dot(n0, n1) < 0)
                //n1 = negate(n1);

            double costheta = fabs(dot(divide(n0,norm(n0)), divide(n1,norm(n1))));

            return std::acos(costheta);
        }

        cart3_point  rotate_around_axis(
                const cart3_polygon& dst, 
                unsigned int a0, unsigned int a1,
                unsigned int pt_idx,
                double angle){

            using namespace boost::geometry::strategy::transform;
            cart3_point rotax = difference(dst.outer()[a1], dst.outer()[a0]);
            rotax = divide(rotax, norm(rotax));

            cart3_point orig = dst.outer()[a0];
            cart3_point ptx = dst.outer()[pt_idx];
            ptx = difference(ptx, orig);
            double pt[3] = {ptx.get<0>(), ptx.get<1>(), ptx.get<2>()};

            double x = rotax.get<0>(), y = rotax.get<1>(), z = rotax.get<2>();
            double c = cos(angle), s = sin(angle);
            double x2 = x*x, y2 = y*y, z2 = z*z;
            double rmat[3][3] = {
                {c+x2*(1-c),    x*y*(1-c)-z*s, x*z*(1-c)+y*s},
                {y*x*(1-c)+z*s, c+y2*(1-c), y*z*(1-c)-x*s},
                {z*x*(1-c)-y*s, z*y*(1-c)+x*s, c+z2*(1-c)}};

            cart3_point pt2(
                    rmat[0][0] * pt[0] + rmat[0][1] * pt[1] + rmat[0][2] * pt[2],
                    rmat[1][0] * pt[0] + rmat[1][1] * pt[1] + rmat[1][2] * pt[2],
                    rmat[2][0] * pt[0] + rmat[2][1] * pt[1] + rmat[2][2] * pt[2]);

            pt2 = sum(pt2, orig);
            
            return pt2;
        }
        namespace detail
        {
            template<class T>
            shared_edge_property determine_shared_vertices(const T& a, const T& b){
                assert(a.outer().size() == 3);
                assert(b.outer().size() == 3);
                shared_edge_property sep;

                int a2b[3] = {0, 0, 0};
                int n_found = 0;
                for (int i = 0; i < 3; ++i)
                {
                    bool found = false;
                    for (int j = 0; j < 3; ++j)
                    {
                        if(norm2(difference(a.outer()[i], b.outer()[j])) < 0.000001)
                        {
                            a2b[i] = j;
                            found = true;
                            n_found ++;
                            break;
                        }
                    }
                    if(!found){
                        sep.single_l = i;  // the point in b which is not in a
                        sep.shared_l0 = (i + 1) % 3;  // another point.
                        sep.shared_l1 = (i + 2) % 3;  // a0 and a1 are the rotation axis
                    }
                }
                if(n_found != 2){
                    BOOST_LOG_TRIVIAL(error) << "n_found: " << n_found;
                }
                    
                assert(n_found == 2);
                sep.shared_h0 = a2b[sep.shared_l0];
                sep.shared_h1 = a2b[sep.shared_l1];
                if(sep.shared_h0 != 0 && sep.shared_h1 != 0)
                    sep.single_h = 0;
                else if(sep.shared_h0 != 1 && sep.shared_h1 != 1)
                    sep.single_h = 1;
                else if(sep.shared_h0 != 2 && sep.shared_h1 != 2)
                    sep.single_h = 2;
                else{
                    assert(false);
                }
                assert(sep.shared_h0 != sep.shared_h1);
                assert(sep.shared_l0 != sep.shared_l1);
                assert(sep.single_h != sep.shared_h1);
                assert(sep.single_l != sep.shared_l1);
                assert(sep.single_h != sep.shared_h0);
                assert(sep.single_l != sep.shared_l0);
                return sep;
            }
        }
        shared_edge_property determine_shared_vertices(const cart3_polygon& a, const cart3_polygon& b){
            return detail::determine_shared_vertices(a, b);
        }
        shared_edge_property determine_shared_vertices(const polar2_polygon& a, const polar2_polygon& b){
            return detail::determine_shared_vertices(a, b);
        }

        cart3_polygon flatten(const cart3_polygon& t0c, const cart3_polygon& t1c,
                const cart3_polygon& base, const shared_edge_property& , double fact){
            using namespace boost::geometry::strategy::transform;

            auto sep = determine_shared_vertices(t0c, t1c);
            cart3_polygon dst = construct_triangle_pair(t0c, t1c, base, sep);

            unsigned char a10 = sep.shared_h0,
                          a11 = sep.shared_h1,
                          pt1 = sep.single_h;

            double ang = angle_between(base, dst);
            cart3_point pt1_orig = dst.outer()[pt1];
            {
                if(fact != 0.){
                    dst.outer()[pt1] = rotate_around_axis(dst, a10, a11, pt1, fact * ang);
                    double ang2 = angle_between(base, dst);
                    if(fabs(ang2) > fabs(ang)){
                        dst.outer()[pt1] = pt1_orig;
                        dst.outer()[pt1] = rotate_around_axis(dst, a10, a11, pt1, -fact * ang);
                    }
                    //double ang3 = angle_between(base, dst);
                    //std::cout << "New angle betwen base and dst: " << ang3 << std::endl;
                }
            }
            return dst;
        }
        cart3_point relative_tria_pos(const cart3_point& d0, const cart3_point& d1, const cart3_point& n0, const cart3_point p){
            namespace ublas = boost::numeric::ublas;
            namespace atlas = boost::numeric::bindings::atlas;
            typedef ublas::matrix<double, ublas::column_major> matrix;
            std::vector<int> ipiv(3);

            // A x = b,   rhs is b [input] and x [output]
            matrix A(3,3);
            matrix rhs(3, 1);

            rhs(0,0) = p.get<0>();
            rhs(1,0) = p.get<1>();
            rhs(2,0) = p.get<2>();

            A(0, 0) = d0.get<0>();
            A(1, 0) = d0.get<1>();
            A(2, 0) = d0.get<2>();
            A(0, 1) = d1.get<0>();
            A(1, 1) = d1.get<1>();
            A(2, 1) = d1.get<2>();
            A(0, 2) = n0.get<0>();
            A(1, 2) = n0.get<1>();
            A(2, 2) = n0.get<2>();

            int res0 = atlas::getrf(A, ipiv);
            int res1 = atlas::getrs(A, ipiv, rhs);
            if(res0 != 0){
                BOOST_LOG_TRIVIAL(warning) << "res0 != 0, A: " << A;
            }
            if(res1 != 0){
                BOOST_LOG_TRIVIAL(warning) << "res1 != 0, A: " << A;
            }
            assert(res0 == 0);
            assert(res1 == 0);
            return cart3_point(rhs(0, 0), rhs(1,0), rhs(2,0));
        }
        cart3_point relative_tria_pos(const cart3_polygon& t, const cart3_point& p){
            cart3_point d0 = difference(t.outer()[1], t.outer()[0]);
            cart3_point d1 = difference(t.outer()[2], t.outer()[0]);
            cart3_point n0 = cross(d0, d1);
            return relative_tria_pos(d0, d1, n0, difference(p, t.outer()[0]));
        }
        cart3_point reltriapos_to_abs(const cart3_point& d0, const cart3_point& d1, const cart3_point& n0, const cart3_point p){
            namespace ublas = boost::numeric::ublas;
            namespace atlas = boost::numeric::bindings::atlas;
            typedef ublas::matrix<double, ublas::column_major> matrix;
            matrix rhs(3, 1);
            matrix A(3, 3);

            rhs(0,0) = p.get<0>();
            rhs(1,0) = p.get<1>();
            rhs(2,0) = p.get<2>();

            A(0, 0) = d0.get<0>();
            A(1, 0) = d0.get<1>();
            A(2, 0) = d0.get<2>();
            A(0, 1) = d1.get<0>();
            A(1, 1) = d1.get<1>();
            A(2, 1) = d1.get<2>();
            A(0, 2) = n0.get<0>();
            A(1, 2) = n0.get<1>();
            A(2, 2) = n0.get<2>();

            rhs = ublas::prod(A, rhs);
            return cart3_point(rhs(0,0), rhs(1,0), rhs(2,0));
        }
        cart3_point reltriapos_to_abs(const cart3_polygon& t, const cart3_point& p){
            cart3_point d0 = difference(t.outer()[1], t.outer()[0]);
            cart3_point d1 = difference(t.outer()[2], t.outer()[0]);
            cart3_point n0 = cross(d0, d1);
            return sum(t.outer()[0], reltriapos_to_abs(d0, d1, n0, p));
        }
        cart3_polygon construct_triangle_pair(
                const cart3_polygon& t0, const cart3_polygon& t1,
                const cart3_polygon& base,
                const shared_edge_property& sep
                ){
            cart3_polygon dst;
            unsigned char a00 = sep.shared_l0,
                          a01 = sep.shared_l1,
                          pt0 = sep.single_l;
            unsigned char a10 = sep.shared_h0,
                          a11 = sep.shared_h1,
                          pt1 = sep.single_h;
            dst.outer().resize(3);
            dst.outer()[a10] = base.outer()[a00];
            dst.outer()[a11] = base.outer()[a01];

            // express pt0 as a linear combination of d00, d01, and n0
            cart3_point d0 = difference(t0.outer()[a00], t0.outer()[pt0]);
            cart3_point d1 = difference(t0.outer()[a01], t0.outer()[pt0]);
            cart3_point d0b = difference(base.outer()[a00], base.outer()[pt0]);
            cart3_point d1b = difference(base.outer()[a01], base.outer()[pt0]);

            cart3_point n0 = cross(
                    difference(t0.outer()[1], t0.outer()[0]), 
                    difference(t0.outer()[2], t0.outer()[0]));

            cart3_point n0b = cross(
                    difference(base.outer()[1], base.outer()[0]), 
                    difference(base.outer()[2], base.outer()[0]));

            namespace ublas = boost::numeric::ublas;
            namespace atlas = boost::numeric::bindings::atlas;
            typedef ublas::matrix<double, ublas::column_major> matrix;
            std::vector<int> ipiv(3);

            // A x = b,   rhs is b [input] and x [output]
            matrix A(3,3);
            matrix rhs(3, 1);

            cart3_point pt1x = difference(t1.outer()[pt1], t0.outer()[pt0]);

            rhs(0,0) = pt1x.get<0>();
            rhs(1,0) = pt1x.get<1>();
            rhs(2,0) = pt1x.get<2>();

            A(0, 0) = d0.get<0>();
            A(1, 0) = d0.get<1>();
            A(2, 0) = d0.get<2>();
            A(0, 1) = d1.get<0>();
            A(1, 1) = d1.get<1>();
            A(2, 1) = d1.get<2>();
            A(0, 2) = n0.get<0>();
            A(1, 2) = n0.get<1>();
            A(2, 2) = n0.get<2>();

            int res0 = atlas::getrf(A, ipiv);
            int res1 = atlas::getrs(A, ipiv, rhs);
            if(res0 != 0){
                std::cout << "res0 != 0, A: " << A << std::endl;
            }
            if(res1 != 0){
                std::cout << "res1 != 0, A: " << A << std::endl;
            }
            assert(res0 == 0);
            assert(res1 == 0);

            A(0, 0) = d0b.get<0>();
            A(1, 0) = d0b.get<1>();
            A(2, 0) = d0b.get<2>();
            A(0, 1) = d1b.get<0>();
            A(1, 1) = d1b.get<1>();
            A(2, 1) = d1b.get<2>();
            A(0, 2) = n0b.get<0>();
            A(1, 2) = n0b.get<1>();
            A(2, 2) = n0b.get<2>();

            rhs = ublas::prod(A, rhs);
            //std::cout << "A1: " << A << std::endl;
            //std::cout << "rhs: "<<rhs(0,0) <<" " << rhs(1,0) << " " << rhs(2,0)<<std::endl;

            // determine pt1 using the same linear combination
            dst.outer()[pt1].set<0>(rhs(0,0));
            dst.outer()[pt1].set<1>(rhs(1,0));
            dst.outer()[pt1].set<2>(rhs(2,0));

            dst.outer()[pt1] = sum(dst.outer()[pt1], base.outer()[pt0]);
            return dst;
        }
        bool inside_box(const polar2_polygon& sp, double lat_off, double lng_off){
            bool inside = true;
            for(const auto& p : sp.outer()){
                if(fabs(p.get<0>()) > lat_off || fabs(p.get<1>()) > lng_off) {
                    inside = false;
                    break;
                }
            }
            return inside;
        }
        polar2_point spherical_centroid(const polar2_polygon& poly){
            using namespace boost::geometry::strategy::transform;
            assert(poly.outer().size() == 3);
            polar2_point res;
            cart3_polygon c2p;
            cart3_point c2pt(0., 0., 0.);
            from_spherical_equatorial_2_to_cartesian_3<polar2_point, cart3_point> sph2cart3;
            from_cartesian_3_to_spherical_equatorial_2<cart3_point, polar2_point> cart3sph2;
            bg::transform(poly, c2p, sph2cart3);
            for(auto p : c2p.outer()){
                c2pt.set<0>( c2pt.get<0>() + p.get<0>() );
                c2pt.set<1>( c2pt.get<1>() + p.get<1>() );
                c2pt.set<2>( c2pt.get<2>() + p.get<2>() );
            }
            c2pt.set<0>( c2pt.get<0>() / 3. );
            c2pt.set<1>( c2pt.get<1>() / 3. );
            c2pt.set<2>( c2pt.get<2>() / 3. );
            bg::transform(c2pt, res, cart3sph2);
            return res;
        }
    cart3_point 
        tria2tria(const cart3_polygon& t0, const cart3_polygon& t1,
                const cart3_point& src){
            return reltriapos_to_abs(t1, relative_tria_pos(t0, src));
        }
    cart3_polygon
        tria2tria(const cart3_polygon& t0, const cart3_polygon& t1,
                const cart3_polygon& src){
            cart3_polygon ret;
            for(auto p : src.outer())
                ret.outer().push_back(reltriapos_to_abs(t1, relative_tria_pos(t0, p)));
            return ret;
        }

    }
}
