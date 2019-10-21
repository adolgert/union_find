/*
 * morton.hpp
 *
 *  Created on: Mar 21, 2012
 *      Author: ajd27
 */

#ifndef MORTON_HPP_
#define MORTON_HPP_

#include <boost/array.hpp>


namespace raster_stats
{

    /*! This creates a size_t with val repeats of the two-bit pattern in base.
     *  It exists to make masks for Morton order that look like
     *  0b010101 and 0b101010, except 64 bits long.
     *  xmask=alternating_bits<32,0b01>::value  => 0101010101...
     *  ymask=alternating_bits<32,0b10>::value  => 1010101010...
     */
    template<int val,size_t base>
    struct alternating_bits {
        typedef size_t value_type;
        static value_type const value=(alternating_bits<val-1,base>::value <<2)+base;
    };

    template<size_t base>
    struct alternating_bits<1,base>
    {
        typedef size_t value_type;
        static value_type const value=base;
    };

    template<size_t a,size_t b>
    struct tmin
    {
        typedef size_t value_type;
        static value_type const value=(a>b?b:a);
    };


    /*! This preprocesses a constant morton value.
     *  The morton_d will treat x as the dth dimension.
     *  d is the dimension of this coordinate x.
     *  D is the total dimension of the space (2 for x and y).
     */
    template<size_t x, size_t d, size_t D>
    struct morton_d
    {
        typedef size_t value_type;
        static value_type const value=(morton_d<(x>>1),d,D>::value <<D)|((x&1)<<d);
    };

    template<size_t d,size_t D>
    struct morton_d<0,d,D>
    {
      typedef size_t value_type;
      static value_type const value=0;
    };

    template<size_t d,size_t D>
    struct morton_d<1,d,D>
    {
      typedef size_t value_type;
      static value_type const value=(1<<d);
    };


    template<size_t x,size_t y>
    struct morton_xy
    {
        typedef size_t value_type;
        static value_type const value=(morton_d<x,0,2>::value|morton_d<y,1,2>::value);
    };



    struct morton_calculations
    {
        /*! Combine an x and y cooordinate into a single Morton coordinate.
         *  The types of x and y will only use half the bits of the Morton coord
         *  so they may be represented with a smaller type.
         */
        template<typename XY, typename M>
        static M combine_xy(const boost::array<XY,2>& x)
        {
            const M relevant_bits = tmin<sizeof(M)/2,sizeof(XY)>::value * CHAR_BIT;
            M z=0;
            for (M i=0; i < relevant_bits; i++) {
                M bit = 1u << i;
                for (size_t dim=0; dim<2; dim++) {
                    z |= (x[dim] & bit) << (i+dim);
                }
            }
            return z;
        }


        template<typename XY, typename M>
        static boost::array<XY,2> detangle(M n)
        {
            const XY relevant_bits = tmin<sizeof(M)/2,sizeof(XY)>::value * CHAR_BIT;
            boost::array<XY,2> xy = {{0, 0}};
            for (XY i=0; i<relevant_bits; i++) {
                for (size_t dim=0; dim<2; dim++) {
                    xy[dim] |= ( n & (1<<(2*i+dim)) ) >> (i+dim);
                }
            }
            return xy;
        }


        /*! Given an x and y interleaved into a single value, add another interleaved value.
         *  A Morton coordinate, added to an interleaved x, will increment only in x.
         *  This is basically a bit-shifting re-implementation of addition over every
         *  other bit.
         */
        template<typename M>
        static M add_interleaved(M n, M dn)
        {
            while (n) {
                M a = n & dn;
                M b = n ^ dn;
                n = a << 2;
                dn = b;
            }
            return dn;
        }


    };
    }

#endif /* MORTON_HPP_ */
