/*
 * Jon Leech (leech @ cs.unc.edu) 3/24/89
 *           http://www.cs.unc.edu/~jon
 * icosahedral code added by Jim Buddenhagen (jb1556@daditz.sbc.com) 5/93
 */

/*% cc -g sphere.c -o sphere -lm
 *
 * sphere - generate a triangle mesh approximating a sphere by
 *  recursive subdivision. First approximation is an platonic
 *  solid; each level of refinement increases the number of
 *  triangles by a factor of 4.
 *
 * Level 3 (128 triangles for an octahedron) is a good tradeoff if
 *  gouraud shading is used to render the database.
 *
 * Usage: sphere [level] [-p] [-c] [-f] [-t] [-i]
 *	level is an integer >= 1 setting the recursion level (default 1).
 *	-p causes generation of a PPHIGS format ASCII archive
 *	    instead of the default generic output format.
 *	-c causes triangles to be generated with vertices in counterclockwise
 *	    order as viewed from the outside in a RHS coordinate system.
 *	    The default is clockwise order.
 *	-f generates triangle without per-vertex normals (PPHIGS only)
 *	-t starts with a tetrahedron instead of an octahedron
 *	-i starts with a icosahedron instead of an octahedron
 *
 *  The subroutines print_object() and print_triangle() should
 *  be changed to generate whatever the desired database format is.
 *

*/
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "myriaworld.h"
#include "SPHEmesh.hpp"
using myriaworld::cart3_point;
using myriaworld::cart3_polygon;

typedef struct {
    double  x, y, z;
} point;

typedef struct {
    double  w0, w1, wc;
} edgeweight;

typedef struct {
    point     pt[3];	/* Vertices of triangle */
    double    area;	/* Unused; might be used for adaptive subdivision */
    edgeweight edges[3]; /* edge properties of triangle */
} triangle;

typedef struct {
    int       npoly;	/* # of triangles in object */
    triangle *poly;	/* Triangles */
} object;

/* Six equidistant points lying on the unit sphere */
#define XPLUS {  1,  0,  0 }	/*  X */
#define XMIN  { -1,  0,  0 }	/* -X */
#define YPLUS {  0,  1,  0 }	/*  Y */
#define YMIN  {  0, -1,  0 }	/* -Y */
#define ZPLUS {  0,  0,  1 }	/*  Z */
#define ZMIN  {  0,  0, -1 }	/* -Z */

#define EW {{2, 1, 2}, {1, 1, 2}, {1, 1, 2}}

/* Vertices of a unit octahedron */
triangle octahedron[] = {
    { { XPLUS, ZPLUS, YPLUS }, 0.0, EW },
    { { YPLUS, ZPLUS, XMIN  }, 0.0, EW },
    { { XMIN , ZPLUS, YMIN  }, 0.0, EW },
    { { YMIN , ZPLUS, XPLUS }, 0.0, EW },
    { { XPLUS, YPLUS, ZMIN  }, 0.0, EW },
    { { YPLUS, XMIN , ZMIN  }, 0.0, EW },
    { { XMIN , YMIN , ZMIN  }, 0.0, EW },
    { { YMIN , XPLUS, ZMIN  }, 0.0, EW }
};

/* A unit octahedron */
object oct = {
    sizeof(octahedron) / sizeof(octahedron[0]),
    &octahedron[0]
};

/* Vertices of a tetrahedron */
#define sqrt_3 0.5773502692
#define PPP {  sqrt_3,	sqrt_3,  sqrt_3 }   /* +X, +Y, +Z */
#define MMP { -sqrt_3, -sqrt_3,  sqrt_3 }   /* -X, -Y, +Z */
#define MPM { -sqrt_3,	sqrt_3, -sqrt_3 }   /* -X, +Y, -Z */
#define PMM {  sqrt_3, -sqrt_3, -sqrt_3 }   /* +X, -Y, -Z */

/* Structure describing a tetrahedron */
triangle tetrahedron[] = {
    {{ PPP, MMP, MPM }, 0.0, EW},
    {{ PPP, PMM, MMP }, -1.0, EW},
    {{ MPM, MMP, PMM }, 0.0, EW},
    {{ PMM, PPP, MPM }, 0.0, EW}
};

object tet = {
    sizeof(tetrahedron) / sizeof(tetrahedron[0]),
    &tetrahedron[0]
};

/* Twelve vertices of icosahedron on unit sphere */
#define tau 0.8506508084      /* t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)  */
#define one 0.5257311121      /* one=1/sqrt(1+t^2) , unit sphere     */
#define ZA {  tau,  one,    0 }
#define ZB { -tau,  one,    0 }
#define ZC { -tau, -one,    0 }
#define ZD {  tau, -one,    0 }
#define YA {  one,   0 ,  tau }
#define YB {  one,   0 , -tau }
#define YC { -one,   0 , -tau }
#define YD { -one,   0 ,  tau }
#define XA {   0 ,  tau,  one }
#define XB {   0 , -tau,  one }
#define XC {   0 , -tau, -one }
#define XD {   0 ,  tau, -one }

/* Structure for unit icosahedron */
triangle icosahedron[] = {
    { { YA, XA, YD }, 0.0, EW },
    { { YA, YD, XB }, 0.0, EW },
    { { YB, YC, XD }, 0.0, EW },
    { { YB, XC, YC }, 0.0, EW },
    { { ZA, YA, ZD }, 0.0, EW },
    { { ZA, ZD, YB }, 0.0, EW },
    { { ZC, YD, ZB }, 0.0, EW },
    { { ZC, ZB, YC }, 0.0, EW },
    { { XA, ZA, XD }, 0.0, EW },
    { { XA, XD, ZB }, 0.0, EW },
    { { XB, XC, ZD }, 0.0, EW },
    { { XB, ZC, XC }, 0.0, EW },
    { { XA, YA, ZA }, 0.0, EW },
    { { XD, ZA, YB }, 0.0, EW },
    { { YA, XB, ZD }, 0.0, EW },
    { { YB, ZD, XC }, 0.0, EW },
    { { YD, XA, ZB }, 0.0, EW },
    { { YC, ZB, XD }, 0.0, EW },
    { { YD, ZC, XB }, 0.0, EW },
    { { YC, XC, ZC }, 0.0, EW }
};

/* A unit icosahedron */
object ico = {
    sizeof(icosahedron) / sizeof(icosahedron[0]),
    &icosahedron[0]
};

int PPHIGSflag = 0; /* Don't generate PPHIGS format output */
int Flatflag = 0;   /* Don't generate per-vertex normals */

/* Forward declarations */
point *normalize(point *p);
point *midpoint(point *a, point *b);
void flip_object(object *obj );
void print_object(object *obj);
void print_triangle(triangle *t);
void pphigs_header(int);
void pphigs_trailer();
void split_edge_weights(const edgeweight& src, edgeweight& dst0, edgeweight& dst1){
    dst0.w0 = src.w0;
    dst1.w0 = src.wc;

    dst0.w1 = src.wc;
    dst1.w1 = src.w1;

    dst0.wc = (src.w0 + src.wc) / 2.;
    dst1.wc = (src.wc + src.w1) / 2.;
}

extern char *malloc(/* unsigned */);

std::vector<cart3_polygon>
triangulate_sphere(int maxlevel)
{
    object *old = &ico,		/* Default is octahedron */
           *novel;
    int     ccwflag = 0,	/* Reverse vertex order if true */
            i,
            level;		/* Current subdivision level */

    if (ccwflag)
        flip_object(old);

    /* Subdivide each starting triangle (maxlevel - 1) times */
    for (level = 1; level < maxlevel; level++) {
        /* Allocate a novel object */
        novel = (object *)malloc(sizeof(object));
        if (novel == NULL) {
            fprintf(stderr, "Out of memory on subdivision level %d\n", level);
            exit(1);
        }
        novel->npoly = old->npoly * 4;

        /* Allocate 4* the number of points in the current approximation */
        novel->poly  = (triangle *)malloc(novel->npoly * sizeof(triangle));
        if (novel->poly == NULL) {
            fprintf(stderr, "Out of memory on subdivision level %d\n", level);
            exit(1);
        }

        /* Subdivide each triangle in the old approximation and normalize
         *  the novel points thus generated to lie on the surface of the unit
         *  sphere.
         * Each input triangle with vertices labelled [0,1,2] as shown
         *  below will be turned into four novel triangles:
         *
         *			Make novel points
         *			    a = (0+2)/2
         *			    b = (0+1)/2
         *			    c = (1+2)/2
         *       1
         *       /\		Normalize a, b, c
         *      /  \
         *    b/____\ c		Construct novel triangles
         *    /\    /\		    [0,b,a]
         *   /  \  /  \		    [b,1,c]
         *  /____\/____\	    [a,b,c]
         * 0      a     2	    [a,c,2]
         */
        for (i = 0; i < old->npoly; i++) {
            triangle
                *oldt = &old->poly[i],
                *newt = &novel->poly[i*4];
            point a, b, c;

            a = *normalize(midpoint(&oldt->pt[0], &oldt->pt[2]));
            b = *normalize(midpoint(&oldt->pt[0], &oldt->pt[1]));
            c = *normalize(midpoint(&oldt->pt[1], &oldt->pt[2]));

            newt[0].pt[0] = oldt->pt[0];
            newt[0].pt[1] = b;
            newt[0].pt[2] = a;

            newt[1].pt[0] = b;
            newt[1].pt[1] = oldt->pt[1];
            newt[1].pt[2] = c;

            newt[2].pt[0] = a;
            newt[2].pt[1] = b;
            newt[2].pt[2] = c;

            newt[3].pt[0] = a;
            newt[3].pt[1] = c;
            newt[3].pt[2] = oldt->pt[2];

            for(unsigned int k=0; k<4; k++){
                for(unsigned int e=0; e<3; e++){
                    newt[k].edges[e].w0 = level;
                    newt[k].edges[e].w1 = level;
                    newt[k].edges[e].wc = level + 1;
                }
            }

            split_edge_weights(oldt->edges[0], newt[0].edges[0], newt[1].edges[0]);
            split_edge_weights(oldt->edges[1], newt[1].edges[1], newt[3].edges[1]);
            split_edge_weights(oldt->edges[2], newt[3].edges[2], newt[0].edges[2]);

        }

        if (level > 1) {
            free(old->poly);
            free(old);
        }

        /* Continue subdividing novel triangles */
        old = novel;
    }

    std::vector<cart3_polygon> ret(old->npoly);
    for (int i = 0; i < old->npoly; i++){
        //print_triangle(&old->poly[i]);
        triangle* t = &old->poly[i];
        for(int j=0; j < 3; j++){
            ret[i].outer().push_back(
                    cart3_point(t->pt[j].x, t->pt[j].y, t->pt[j].z));
        }
    }
    return ret;
}

/* Normalize a point p */
point *normalize(point* p)
{
    static point r;
    double mag;

    r = *p;
    mag = r.x * r.x + r.y * r.y + r.z * r.z;
    if (mag != 0.0) {
        mag = 1.0 / sqrt(mag);
        r.x *= mag;
        r.y *= mag;
        r.z *= mag;
    }

    return &r;
}

/* Return the midpoint on the line between two points */
point *midpoint(point* a, point* b)
{
    static point r;

    r.x = (a->x + b->x) * 0.5;
    r.y = (a->y + b->y) * 0.5;
    r.z = (a->z + b->z) * 0.5;

    return &r;
}

/* Reverse order of points in each triangle */
void flip_object(object* obj)
{
    int i;
    for (i = 0; i < obj->npoly; i++) {
        point tmp;
        tmp = obj->poly[i].pt[0];
        obj->poly[i].pt[0] = obj->poly[i].pt[2];
        obj->poly[i].pt[2] = tmp;
    }
}

/* Write out all triangles in an object */
void print_object(object* obj)
{
    for (int i = 0; i < obj->npoly; i++)
        print_triangle(&obj->poly[i]);
}

/* Output a triangle */
void print_triangle(triangle* t)
{
    for(int i=0; i < 3; i++){
        if(i>0)
            printf("\t");
        printf("\t%3.9f %3.9f %3.9f",
                t->pt[i].x, t->pt[i].y, t->pt[i].z);
    }
    printf("\n");
}
