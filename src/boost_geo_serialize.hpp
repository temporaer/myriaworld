#ifndef __BOOST_GEO_SERIALIZE_HPP__
#     define __BOOST_GEO_SERIALIZE_HPP__


// boost
#include "myriaworld.h"
#include <boost/serialization/split_free.hpp>
//#include <boost/unordered_set.hpp>


namespace boost
{
    /////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8/////////9/////////A
    namespace serialization
    {
        using myriaworld::polar2_point;
        using myriaworld::cart2_point;
        using myriaworld::cart3_point;
        using myriaworld::polar2_polygon;
        using myriaworld::cart2_polygon;
        using myriaworld::cart3_polygon;

        template <class Archive>
            void save(Archive &ar, const myriaworld::polar2_point &pos, const
                    unsigned int )
            {
                double p0(pos.get<0>()), p1(pos.get<1>());
                ar & BOOST_SERIALIZATION_NVP(p0);
                ar & BOOST_SERIALIZATION_NVP(p1);
            }

        template <class Archive>
            void load(Archive &ar, myriaworld::polar2_point &pos, const unsigned
                    int )
            {
                double p0, p1;
                ar & BOOST_SERIALIZATION_NVP(p0);
                ar & BOOST_SERIALIZATION_NVP(p1);
                pos = myriaworld::polar2_point(p0, p1);
            }


        template<class Archive>
            inline void serialize(Archive & ar, boost::geometry::model::ring<polar2_point, false, false> &ring, const unsigned int )
            {
                ar & static_cast<std::vector<polar2_point>& >(ring);
            }
        template<class Archive>
            inline void serialize(Archive & ar, boost::geometry::model::ring<cart2_point, false, false> &ring, const unsigned int )
            {
                ar & static_cast<std::vector<cart2_point>& >(ring);
            }
        template<class Archive>
            inline void serialize(Archive & ar, boost::geometry::model::ring<cart3_point, false, false> &ring, const unsigned int )
            {
                ar & static_cast<std::vector<cart3_point>& >(ring);
            }

                template<class Archive>
                inline void serialize(Archive & ar, polar2_polygon &t, const unsigned int )
                {
                        ar & t.outer();
                        ar & t.inners();
                }
                template<class Archive>
                inline void serialize(Archive & ar, cart2_polygon &t, const unsigned int )
                {
                        ar & t.outer();
                        ar & t.inners();
                }
                template<class Archive>
                inline void serialize(Archive & ar, cart3_polygon &t, const unsigned int )
                {
                        ar & t.outer();
                        ar & t.inners();
                }

        template <class Archive>
            void save(Archive &ar, const myriaworld::cart2_point &pos, const
                    unsigned int )
            {
                double lat(pos.get<0>()), lon(pos.get<1>());
                ar & BOOST_SERIALIZATION_NVP(lon);
                ar & BOOST_SERIALIZATION_NVP(lat);
            }

        template <class Archive>
            void load(Archive &ar, myriaworld::cart2_point &pos, const unsigned
                    int )
            {
                double lat, lon;
                ar & BOOST_SERIALIZATION_NVP(lon);
                ar & BOOST_SERIALIZATION_NVP(lat);
                pos = myriaworld::cart2_point(lon, lat);
            }
        template <class Archive>
            void save(Archive &ar, const myriaworld::cart3_point &pos, const
                    unsigned int )
            {
                double lat(pos.get<0>()), lon(pos.get<1>()), d(pos.get<2>());
                ar & BOOST_SERIALIZATION_NVP(lon);
                ar & BOOST_SERIALIZATION_NVP(lat);
                ar & BOOST_SERIALIZATION_NVP(d);
            }

        template <class Archive>
            void load(Archive &ar, myriaworld::cart3_point &pos, const unsigned
                    int )
            {
                double lat, lon, d;
                ar & BOOST_SERIALIZATION_NVP(lon);
                ar & BOOST_SERIALIZATION_NVP(lat);
                ar & BOOST_SERIALIZATION_NVP(d);
                pos = myriaworld::cart3_point(lon, lat, d);
            }

    } // namespace serialization
    /////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8/////////9/////////A
    namespace geometry
    {
        size_t hash_value(const myriaworld::polar2_point &pnt)
        {
            return boost::hash_value(pnt.get<0>()) +
                boost::hash_value(pnt.get<1>());
        }

        bool operator==(const myriaworld::polar2_point &lhs, const myriaworld::polar2_point&rhs)
        {
            return lhs.get<0>() == rhs.get<0>() && lhs.get<1>() == rhs.get<1>();
        }
    } // namespace geometry

} // namespace boost
BOOST_SERIALIZATION_SPLIT_FREE(myriaworld::polar2_point)
BOOST_SERIALIZATION_SPLIT_FREE(myriaworld::cart2_point)
BOOST_SERIALIZATION_SPLIT_FREE(myriaworld::cart3_point)

#endif /* __BOOST_GEO_SERIALIZE_HPP__ */
