#include <complex>
#include <random>
#include <cmath>
#include <chrono>
#include <iostream>

using C = std::complex<float>;

int main(void)
{
    std::vector<C> input(1024*1024);

    std::random_device rd;
    std::mt19937 e2(rd());
    std::uniform_real_distribution<> dist(0, 1);

    for(int i = 0; i < 30; ++i)
    for(auto& e: input)
        e = C(dist(e2), dist(e2));

    auto start = std::chrono::high_resolution_clock::now();
    for(auto& e: input)
        e = e + std::conj(e) * C(0.1);

    for(auto& e: input){
        auto x = e;
        auto r = std::norm(x);

        r = 1.0 / (sqrt(1 + r));
        e = r * x;
    }

    for(auto& e: input){
        auto x = e;
        auto r = std::norm(x);

        r = 1.0 / (cbrt(sqrt(1 + pow(r, 3))));
        e = r * x;
    }

    for(auto& e: input)
        e = e + std::conj(e) * C(0.1);

    C cpx = C(0, 0);
    for(auto e: input)
        cpx += e;
    auto finish = std::chrono::high_resolution_clock::now();

    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);

    std::cout << cpx << std::endl;
    std::cout << microseconds.count() / 1000 << std::endl;

    return 0;
}
