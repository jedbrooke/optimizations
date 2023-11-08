#include <iostream>
#include <chrono>



const size_t REPS = 1e8;
const size_t N = 16;
void standard_implementation(u_int64_t a[N], u_int64_t b[N]) {
    for (int i = 0; i < N; i++) {
        for (size_t r = 0; r < REPS; r++) {
            a[i] = a[i] << 32;
            a[i] = a[i] / b[i];
        }
    }
}

void alternate_implementation(u_int64_t a[N], u_int64_t b[N]) {
    for (size_t r = 0; r < REPS; r++) {
        for (size_t i = 0; i < N; i++) {
            a[i] = a[i] << 32;
            a[i] = a[i] / b[i];
        }
    }
}

void standard_implementation(double a[N], double b[N]) {
    for (int i = 0; i < N; i++) {
        for (size_t r = 0; r < REPS; r++) {
            a[i] = a[i] * (1 << 30);
            a[i] = a[i] / b[i];
        }
    }
}

void alternate_implementation(double a[N], double b[N]) {
    for (size_t r = 0; r < REPS; r++) {
        for (size_t i = 0; i < N; i++) {
            a[i] = a[i] * (1 << 30) ;
            a[i] = a[i] / b[i];
        }
    }
}





#ifdef __x86_64__
#include "immintrin.h"
#elif __ARM64_ARCH_8__
#include <arm_neon.h>
    void simd_implementation(double a[N], double b[N]) {
        const double c[] = {1 << 30, 1 << 30};
        const float64x2_t vc = vld1q_f64(c);
        float64x2_t va[N/2];
        float64x2_t vb[N/2];

        for (size_t i = 0; i < N; i += 2) {
            va[i/2] = vld1q_f64(&a[i]);
            vb[i/2] = vld1q_f64(&b[i]);
        }

        for (size_t r = 0; r < REPS; r++) {
            for (int i = 0; i < N/2; i++) {
                va[i] = vmulq_f64(va[i], vc);
                va[i] = vdivq_f64(va[i], vb[i]);
            }
        }

        for (size_t i = 0; i < N; i += 2) {
            vst1q_f64(&a[i], va[i/2]);
        }
    }
#elif __OMP_H
#include <omp.h>
void omp_implementation(double a[N], double b[N]) {
    for (size_t r = 0; r < REPS; r++) {
        for (size_t i = 0; i < N; i++) {
            a[i] = a[i] * (1 << 30) ;
            a[i] = a[i] / b[i];
        }
    }
}
#endif

typedef double test;

int main(int, char**) {
    std::srand(std::time(nullptr));

    test a[N];
    test a_copy[N];
    test b[N];
    test b_copy[N];

    for (size_t i = 0; i < N; i++) {
        test r_a = std::rand();
        a[i] = r_a;
        a_copy[i] = r_a;

        test r_b = std::rand();
        b[i] = r_b;
        b_copy[i] = r_b;
    }


    auto t1 = std::chrono::high_resolution_clock::now();

    alternate_implementation(a, b);

    auto t2 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> ms_double = t2 - t1;

    double standard_time = ms_double.count();

    std::cout << "standard time: " << standard_time << " ms" << std::endl;

    std::cout << "-------------------" << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    simd_implementation(a_copy, b_copy);
    t2 = std::chrono::high_resolution_clock::now();

    ms_double = t2 - t1;

    double alternate_time = ms_double.count();

    std::cout << "alternate time: " << ms_double.count() << " ms" << std::endl;

    std::cout << "-------------------" << std::endl;

    std::cout << "speedup: " << (standard_time / alternate_time) << " x" << std::endl;

    // check accuracy
    bool correct = true;
    for (size_t i = 0; i < N; i++) {
//        std::cout << "a[" << i << "]: " << a[i] << std::endl;
        if (a[i] != a_copy[i]) {
            correct = false;
            std::cout << "a[" << i << "] mismatch!" << "\t" << a[i] << " != " << a_copy[i] << std::endl;
        }
    }
    if (correct) {
        std::cout << "results are correct" << std::endl;
    }




}
