#ifndef __GEO_HPP__
#     define __GEO_HPP__
#include "myriaworld.h"

namespace myriaworld { namespace geo {

    /// shorthand for haversine distance function
    inline double haversine_distance(const polar2_point& a, const polar2_point& b){
        return bg::distance(a, b, bg::strategy::distance::haversine<double>());
    }

    /// check whether the triangle crosses the -180/180 lat line
    bool crosses_dateline(const polar2_polygon& poly);

    //////////////////////////////////////////////////////////
    //   Casting polar to cart2d and vice versa
    //////////////////////////////////////////////////////////

    void cast_polar2_equatorial_2_cart_2(cart2_point& dst, const polar2_point& src, double addlat, double addlng=0.);
    void cast_polar2_equatorial_2_cart_2(cart2_polygon& dst, const polar2_polygon& src, double addlat, double addlng=0.);
    void cast_polar2_equatorial_2_cart_2(cart2_multipolygon& dst, const polar2_multipolygon& src, double addlat, double addlng=0.);
    void cast_cart_2_polar2_equatorial_2(polar2_point& dst, const cart2_point& src, double addlat, double addlng=0.);
    void cast_cart_2_polar2_equatorial_2(polar2_polygon& dst, const cart2_polygon& src, double addlat, double addlng=0.);
    void cast_cart_2_polar2_equatorial_2(polar2_multipolygon& dst, const cart2_multipolygon& src, double addlat, double addlng=0.);
    polar2_point rotated(const polar2_point& sp, double lat_off, double lng_off);
    polar2_polygon rotated(const polar2_polygon& sp, double lat_off, double lng_off);

    double dot(const cart3_point& a, const cart3_point& b);
    cart3_point difference(const cart3_point& a, const cart3_point& b);
    cart3_point sum(const cart3_point& a, const cart3_point& b);
    cart3_point cross(const cart3_point& a, const cart3_point& b);
    cart3_point negate(const cart3_point& a);
    cart3_point divide(const cart3_point& a, double div);
    double norm2(const cart3_point& a);
    double norm(const cart3_point& a);
    double  angle_between(
            const cart3_polygon& t0, 
            const cart3_polygon& t1);

    cart3_point  rotate_around_axis(
            const cart3_polygon& dst, 
            unsigned int a0, unsigned int a1,
            unsigned int pt_idx,
            double angle);

    /// given two triangles t0 and t1 which share a common edge described in sep,
    /// and a triangle base representing t0 in another position, construct a
    /// new triangle dst such that the angles and the relative distances are
    /// equal in the pairs (t0,t1) and (base, dst), with the angle adjusted by fact.
    cart3_polygon 
        flatten(const cart3_polygon& t0, const cart3_polygon& t1,
                const cart3_polygon& base, const shared_edge_property& sep, double fact);

    /// given two triangles t0 and t1 which share a common edge described in sep,
    /// and a triangle base representing t0 in another position, construct a
    /// new triangle dst such that the angles and the relative distances are
    /// equal in the pairs (t0,t1) and (base, dst).
    cart3_polygon construct_triangle_pair(
            const cart3_polygon& t0, const cart3_polygon& t1,
            const cart3_polygon& base, const shared_edge_property& sep);

    bool inside_box(const polar2_polygon& sp, double lat_off, double lng_off);
    polar2_point spherical_centroid(const polar2_polygon& poly);

    } }
#endif /* __GEO_HPP__ */
