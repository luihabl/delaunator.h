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

    std::cout << "finished" << std::endl;

    return 0;
}