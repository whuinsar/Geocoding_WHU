#pragma once

#include <vector>
#include <cassert>
#include <tuple>
#include <functional>
#include <iostream>
#include <array>

#include <Eigen/Dense>
#include <fmt/format.h>

// #define DEBUG

namespace WHU {
	template <size_t N, typename U = double>
	struct PolynomialNT {
		using interpolateResultType = std::array<double, N>;

		Eigen::MatrixXd parameters; // from a_0 to a_{n-1}, where n is the number of degrees.
		double offset_t, zoom_t; // real_t -> [-1, 1] : (real_t  - offset_t) / zoom_t.
		Eigen::Matrix<double, N, 2> transform_xyz; // (real_{xyz} - offset_{xyz}) / zoom_{xyz} \in [-1, 1]
		Eigen::MatrixXd polynomial_coefficients;
		size_t degree;
		U center; // To record year/month/day

		auto at(double t) const {
			interpolateResultType results;
			results.fill(0);

			t = (t - offset_t) / zoom_t;
			double base = 1;
			for (size_t i = 0; i <= degree; ++i) {
				for (size_t j = 0; j < N; ++j) {
					results[j] += base * polynomial_coefficients(j, i);
				};
				base *= t;
			};
			for (size_t i = 0; i < N; ++i) {
				results[i] = results[i] * transform_xyz(i, 1) + transform_xyz(i, 0);
			}

			return results;
		};

		double distance(double other) {
			return std::abs(other - center);
		};

		// degree_of_polynomial = count - 1
		template <typename T>
		static PolynomialNT fromObjects(const std::vector<T>& objects, size_t start, size_t count) {
			assert(start + count <= objects.size());

			PolynomialNT p;
			p.degree = count - 1;
			// TODO: give another mapping rule, to map x,y,z coordinate to [-1, 1]
			auto getOffsetZoom = [&](std::function<double(const T&)> accessor) -> std::tuple<double, double> {
				double sum = 0, min_value = std::numeric_limits<double>::infinity(), max_value = -std::numeric_limits<double>::infinity();

				for (size_t i = start; i < start + count; ++i) {
					sum += accessor(objects[i]);
					min_value = std::min(min_value, accessor(objects[i]));
					max_value = std::max(max_value, accessor(objects[i]));
				}

				return { (max_value + min_value) / 2, (max_value - min_value) / 2 };
			};
			auto get_t = T::getTimeAccessor();
			std::tie(p.offset_t, p.zoom_t) = getOffsetZoom(get_t);

			auto accessors = T::getDataAccessors();
			for (size_t i = 0; i < N; ++i) {
				std::tie(p.transform_xyz(i, 0), p.transform_xyz(i, 1)) = getOffsetZoom(accessors[i]);
			}
			// std::cout << p.transform_xyz << std::endl;

			p.polynomial_coefficients = Eigen::MatrixXd(N, count);
			Eigen::MatrixXd tmp(count, count);
			for (size_t i = start; i < start + count; ++i) {
				double base = 1, multi = (get_t(objects[i]) - p.offset_t) / p.zoom_t;
				for (size_t j = 0; j < count; ++j) {
					tmp(i - start, j) = base;
					base *= multi;
				};
			};

			tmp = tmp.inverse();
			for (size_t i = 0; i < N; ++i) {
				Eigen::MatrixXd b(count, 1);
				for (size_t j = start; j < count + start; ++j) {
					b(j - start, 0) = (accessors[i](objects[j]) - p.transform_xyz(i, 0)) / p.transform_xyz(i, 1);
				};
				p.polynomial_coefficients.block(i, 0, 1, count) = (tmp * b).transpose();

				// std::cout << b << "\n\n";
				// std::cout << tmp.inverse() * p.polynomial_coefficients.block(i, 0, 1, count).transpose() - b;
			};

			if (std::isnan(p.polynomial_coefficients.sum())) {
				/*
				std::cout << p.transform_xyz << "\n\n";
				std::cout << p.polynomial_coefficients << "\n\n";
				std::cout << tmp << "\n";
				*/
				throw std::runtime_error("Failed to export valid coefficient polynomials.\n");
			}
#ifdef DEBUG
			for (size_t i = start; i < start + count; ++i) {
				auto tmp = p.interp(get_t(objects[i]));
				fmt::print("Pts {}.\n           src.                recompute.               diff.\n", i);
				for (size_t j = 0; j < N; ++j) {
					fmt::print("{: >20.8f}\t{: >20.8f}\t{: >20.8f}.\n", accessors[j](objects[i]), tmp[j], accessors[j](objects[i]) - tmp[j]);
				};
			};
#endif 

			p.center = get_t(objects[start]) + (get_t(objects[start + count - 1]) - get_t(objects[start])) / 2;
			return p;
		};

	};

	template<typename T>
	struct MultiPolynomial {
		std::vector<T> polynomial_array;

		template<typename U>
		static MultiPolynomial fromObjects(const std::vector<U>& objects, size_t count) {
			MultiPolynomial tmp;
			for (size_t i = 0; ; ++i) {
				if (i + count > objects.size()) {
					break;
				}

				try {
					tmp.polynomial_array.push_back(T::fromObjects(objects, i, count));
				}
				catch (...) {
					fmt::print("Failed to establish valid polynomial at {}.\n", i);
				}
			}

			return tmp;
		};

		template<typename U>
		auto at(const U& t) {
			size_t index = 0;
			double diff = std::abs(t - polynomial_array[0].center);
			for (size_t i = 1; i < polynomial_array.size(); ++i) {
				double tmp = std::abs(t - polynomial_array[i].center);
				if (tmp < diff) {
					index = i;
					diff = tmp;
				};
			};

			return polynomial_array[index].at(t);
		}
	};
};
