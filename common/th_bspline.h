#include <vector>

// TODO: CRITICAL! These functions are directly dumped from GAMMA.

namespace WHU {
	struct BSpline2D {
		double** coeffs;
		size_t width, height;

		double at(double i, double j) const;

		BSpline2D(const std::vector<double>& source, size_t width);

		~BSpline2D() {
			free(coeffs[0]);
			free(coeffs);
		};
	};
}