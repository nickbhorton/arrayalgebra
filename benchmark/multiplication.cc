#include <cassert>
#include <chrono>
#include <cstdlib>

#include "arrayalgebra.h"

static auto frand() -> double { return static_cast<double>(rand()) / RAND_MAX; }

int main()
{
    aa::mat4 const t1{
        {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}, {13, 14, 15, 16}}
    };
    aa::mat4 r1 = t1 * t1;
    assert(r1[0][0] == 90);
    assert(r1[0][1] == 100);
    assert(r1[0][2] == 110);
    assert(r1[0][3] == 120);

    assert(r1[1][0] == 202);
    assert(r1[1][1] == 228);
    assert(r1[1][2] == 254);
    assert(r1[1][3] == 280);

    assert(r1[2][0] == 314);
    assert(r1[2][1] == 356);
    assert(r1[2][2] == 398);
    assert(r1[2][3] == 440);

    assert(r1[3][0] == 426);
    assert(r1[3][1] == 484);
    assert(r1[3][2] == 542);
    assert(r1[3][3] == 600);

    size_t constexpr thousands = 1000;
    long time{};
    for (size_t i = 0; i < 1000 * thousands; i++) {
        aa::mat4 const m1{
            {{frand(), frand(), frand(), frand()},
             {frand(), frand(), frand(), frand()},
             {frand(), frand(), frand(), frand()},
             {frand(), frand(), frand(), frand()}}
        };
        aa::mat4 const m2{
            {{frand(), frand(), frand(), frand()},
             {frand(), frand(), frand(), frand()},
             {frand(), frand(), frand(), frand()},
             {frand(), frand(), frand(), frand()}}
        };
        auto start = std::chrono::high_resolution_clock::now();
        aa::mat4 const result = m1 * m2;
        auto stop = std::chrono::high_resolution_clock::now();
        time +=
            std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start)
                .count();

        assert(result[0] == result[0]);
    }
    std::cout << static_cast<double>(time) / (1000 * thousands) << "\n";
}
