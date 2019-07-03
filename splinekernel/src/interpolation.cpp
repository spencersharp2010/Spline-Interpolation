#include "interpolation.hpp"
#include "basisfunctions.hpp"
#include "linalg.hpp"

#include <cmath>
#include <exception>

namespace cie
{
	namespace splinekernel
	{

		ControlPointsAndKnotVector interpolateWithBSplineCurve(const ControlPoints2D& interpolationPoints,
			size_t polynomialDegree)
		{
			// Throw exception if number of x-values is not equal number of y-values
			if (interpolationPoints[1].size() != interpolationPoints[0].size())
			{
				throw std::runtime_error("error");
			}

			size_t numberofInterpolationPoints = interpolationPoints[0].size();
			std::vector<double> t_bar(numberofInterpolationPoints);
			t_bar = centripetalParameterPositions(interpolationPoints);

			std::vector<double> knotVector(numberofInterpolationPoints + polynomialDegree + 1);
			knotVector = knotVectorUsingAveraging(t_bar, polynomialDegree);

			linalg::Matrix A(numberofInterpolationPoints, numberofInterpolationPoints, 0);
			ControlPoints2D controlPoints;
			controlPoints[0].reserve(numberofInterpolationPoints);
			controlPoints[1].reserve(numberofInterpolationPoints);

			A(0, 0) = 1;
			A(numberofInterpolationPoints - 1, numberofInterpolationPoints - 1) = 1;

			for (size_t i = 1; i < numberofInterpolationPoints - 1; i++)
			{
				for (size_t j = 0; j < numberofInterpolationPoints; j++)
				{
					double N = evaluateBSplineBasis(t_bar[i], j, polynomialDegree, knotVector);
					A(i, j) = N;
				}

			}
			controlPoints[0] = linalg::solve(A, interpolationPoints[0]);
			controlPoints[1] = linalg::solve(A, interpolationPoints[1]);

			return { {controlPoints }, {knotVector } };
		}

		std::vector<double> centripetalParameterPositions(const ControlPoints2D& interpolationPoints)
		{
			size_t numberOfInterpolationPoints = interpolationPoints[0].size();

			//initialize vector of doubles to store distance between each interpolation point
			std::vector<double> distance(numberOfInterpolationPoints - 1);

			//initialize variable to store summation of all distances
			double totalDistance = 0;

			//loop over number of interpolation points minus one to calculate distance between each
			//interpolation point and also total distance
			for (size_t i = 0; i < numberOfInterpolationPoints - 1; i++)
			{
				double differenceX = interpolationPoints[0][i + 1] - interpolationPoints[0][i];
				double differenceY = interpolationPoints[1][i + 1] - interpolationPoints[1][i];

				distance[i] = sqrt(sqrt(pow(differenceX, 2) + pow(differenceY, 2)));
				totalDistance += distance[i];

			}

			//initialize vector of doubles to store parameter positions
			std::vector<double> t_bar(numberOfInterpolationPoints);

			//initialize first vector element to 0
			t_bar[0] = 0;

			//loop over all interpolation points - 1 to calculate parameter positions
			for (size_t i = 1; i < numberOfInterpolationPoints; i++)
			{
				t_bar[i] = t_bar[i - 1] + distance[i - 1] / totalDistance;
			}

			return { t_bar };
		}

		std::vector<double> knotVectorUsingAveraging(const std::vector<double>& parameterPositions,
			size_t polynomialDegree)
		{
			double p = polynomialDegree;
			size_t n = parameterPositions.size();
			size_t m = n + polynomialDegree + 1;


			//try
		   // {
				// Throw exception if polynomial degree is too high for given number of points
			if ((p + 1) * 2 > m)
			{
				throw std::runtime_error("error");
			}
			// }
			// catch (...)
			 //{
			 //    std::cout << "Polynomial degree is too large." << std::endl;
			 //}

			std::vector<double> knotVector(m);

			for (size_t i = 0; i < polynomialDegree + 1; i++)
			{
				knotVector[i] = 0.0;
			}

			for (size_t i = polynomialDegree + 1; i < n; i++)
			{
				for (size_t j = 1; j < polynomialDegree + 1; j++)
				{
					knotVector[i] += (1 / p) * parameterPositions[i - p - 1 + j];
				}
			}

			for (size_t i = n; i < m; i++)
			{
				knotVector[i] = 1.0;
			}


			return { knotVector };


		}

	} // namespace splinekernel
} // namespace cie