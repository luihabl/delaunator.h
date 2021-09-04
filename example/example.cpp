#include <delaunator.h>

#include <iostream>


int main() {   
    std::vector<float> points = {
        -12.0f, 8.0f,
        -10.0f, 6.0f,
        -12.0f, 4.0f,
         12.0f, 8.0f,
         10.0f, 6.0f,
         12.0f, 4.0f
    };
    
    Delaunator d(points);

    std::cout << "Triangles:" << std::endl;
    for(std::size_t i = 0; i < d.trianglesLen; i+=3) {
        std::cout << d.triangles[i] << "\t" << d.triangles[i+1] << "\t" << d.triangles[i+2] << std::endl;
    }

    return 0;
}