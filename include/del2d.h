#pragma once

#include <array>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

namespace del2d
{

// Type definitions

typedef float fp; // You can change the floating-point precision here
static constexpr std::size_t NIL = (std::numeric_limits<std::size_t>::max)();

// Utilities

inline fp distance(fp x0, fp y0, fp x1, fp y1)
{
    return sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0));
}

inline fp sqr_distance(fp x0, fp y0, fp x1, fp y1)
{
    return (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0);
}


class Mesh
{

public:

    size_t n = 0; // number of points
    std::vector<fp> coords;

    std::vector<size_t> triangles;
    std::vector<size_t> half_edges;

    std::vector<size_t> hull_prev;
    std::vector<size_t> hull_next;
    std::vector<size_t> hull_tri;

    size_t hash_size;
    std::vector<size_t> hull_hash;

    std::vector<size_t> ids;
    std::vector<fp> dists;


    Mesh(const std::vector<fp>& _coords) : coords(_coords) 
    {
        n = coords.size() / 2;
        
        // arrays that will store the triangulation graph
        const size_t max_triangles = std::max(2 * n - 5, 0lu);
        triangles.resize(max_triangles * 3); // is resize or reserve better here?
        half_edges.resize(max_triangles * 3);

        // temporary arrays for tracking the edges of the advancing convex hull
        hash_size = static_cast<std::size_t>(std::ceil(std::sqrt(n)));
        hull_prev.resize(n);
        hull_next.resize(n);
        hull_tri.resize(n);
        
        hull_hash.resize(hash_size); 
        std::fill(hull_hash.begin(), hull_hash.end(), NIL);

        // temporary arrays for sorting points
        ids.resize(n);
        dists.resize(n);
    }

    inline void triangulate()
    {

        // populate an array of point indices; calculate input data bbox
        fp max_x = std::numeric_limits<fp>::lowest();
        fp max_y = std::numeric_limits<fp>::lowest();
        fp min_x = std::numeric_limits<fp>::max();
        fp min_y = std::numeric_limits<fp>::max();

        for (size_t i = 0; i < n; i++)
        {
            const auto x = coords[i * 2 + 0];
            const auto y = coords[i * 2 + 1];

            min_x = std::min(x, min_x);
            min_y = std::min(y, min_y);
            max_x = std::max(x, max_x);
            max_y = std::max(y, max_y);
        
            ids[i] = i;
        }

        const fp span = (max_x - min_x) * (max_x - min_x) + (max_y - min_y) * (max_y - min_y);
        const fp cx = (min_x + max_x) / 2;
        const fp cy = (min_y + max_y) / 2;

        // pick a seed point close to the center
        fp min_dist = std::numeric_limits<fp>::max();
        size_t i0, i1, i2;

        // for(size_t i = 0; i < n; i++)
        

    }





};




















}

