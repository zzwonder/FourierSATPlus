#include <iostream>
#include <chrono>

void function()
{
    long long number = 0;

    for( long long i = 0; i != 2000000000; ++i )
    {
       number += 5;
    }
}

int main()
{
    auto t1 = std::chrono::high_resolution_clock::now();
    function();
    auto t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();

    std::cout << duration<<std::endl;
    return 0;
}
