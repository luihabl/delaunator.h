#include <del2d.h>
#include <iostream>

int main()
{   

    std::vector<float> points = {
        -12, 8,
        -10, 6,
        -12, 4,
        12, 8,
        10, 6,
        12, 4
    };

    auto m = del2d::Delaunator(points);

    for(size_t i = 0; i < m.trianglesLen; i+=3)
    {
        std::cout << m.triangles[i] << " "  << m.triangles[i+1] << " " << m.triangles[i+2] << std::endl;
    }

    std::cout << "finished" << std::endl;

    return 0;
}