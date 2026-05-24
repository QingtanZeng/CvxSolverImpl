#include <cstdint>
#include <cmath>

extern "C" {

    double compute_hypot(double a, double b){
        return std::hypot(a,b);
    }

    void scale_array(double* arr, int32_t size, double factor){
    for (int32_t i = 0; i < size; ++i) {
            arr[i] *= factor;
        }
    }
}