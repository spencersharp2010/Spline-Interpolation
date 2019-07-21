#include "interpolation.hpp"
#include "basisfunctions.hpp"
#include "linalg.hpp"

#include <cmath>
#include <exception>

namespace cie
{
	namespace splinekernel
	{
		// Returns the control points for a b-spline curve with given degree that interpolates the given points.
		ControlPointsAndKnotVector interpolateWithBSplineCurve(const ControlPoints2D& interpolationPoints,
															   size_t polynomialDegree)
		{
			// Throw exception if number of x-values is not equal number of y-values
			if (interpolationPoints[1].size() != interpolationPoints[0].size())
			{
				throw std::runtime_error("error");
			}

			// determine number of given interpolation points
			size_t numberofInterpolationPoints = interpolationPoints[0].size();

			// declare vector t_bar with size equal to number of interpolation points
			std::vector<double> t_bar(numberofInterpolationPoints);

			// calculate t_bar by calling function centripetalParameterPositions
			t_bar = centripetalParameterPositions(interpolationPoints);

			// declare knotVector vector of size n + p + 1
			std::vector<double> knotVector(numberofInterpolationPoints + polynomialDegree + 1);

			// calculate knot vector by calling function knotVectorUsingAveraging
			knotVector = knotVectorUsingAveraging(t_bar, polynomialDegree);

			// declare square matrix "A" of size n x n
			linalg::Matrix A(numberofInterpolationPoints, numberofInterpolationPoints, 0);

			// control points is an array consisting of 2 vectors, each of size n
			ControlPoints2D controlPoints;
			controlPoints[0].reserve(numberofInterpolationPoints);
			controlPoints[1].reserve(numberofInterpolationPoints);

			// set top left and bottom right matrix entry = 1
			A(0, 0) = 1;
			A(numberofInterpolationPoints - 1, numberofInterpolationPoints - 1) = 1;

			// loop across rows of matrix A. First and last rows are known
			for (size_t i = 1; i < numberofInterpolationPoints - 1; i++)
			{
				// loop across columns of matrix A
				for (size_t j = 0; j < numberofInterpolationPoints; j++)
				{
					// calculate each entry, N, of matrix A by calling evaluateBSplineBasis function
					// to calculate N, the parameter position (t_bar), column #, p, and knot vector are needed
					double N = evaluateBSplineBasis(t_bar[i], j, polynomialDegree, knotVector);
					
					// assign N to current entry of matrix A
					A(i, j) = N;
				}

			}
			// solve system of equations for the x-components of the control points
			controlPoints[0] = linalg::solve(A, interpolationPoints[0]);

			// solve system of equations for the y-components of the control points
			controlPoints[1] = linalg::solve(A, interpolationPoints[1]);

			//return control points and the knotVector
			return { {controlPoints }, {knotVector } };
		}

		// function to calculate t_bar vector using centripetal technique
		std::vector<double> centripetalParameterPositions(const ControlPoints2D& interpolationPoints)
		{
			// calculate n by finding size of user-provided interpolationPoints vector
			size_t numberOfInterpolationPoints = interpolationPoints[0].size();

			//initialize vector of doubles to store distance between each interpolation point
			std::vector<double> distance(numberOfInterpolationPoints - 1);

			//initialize variable to store summation of all distances
			double totalDistance = 0;

			//initialize vector of doubles to store parameter positions
			std::vector<double> t_bar(numberOfInterpolationPoints);

			//initialize first vector element to 0
			t_bar[0] = 0;

			//loop over number of interpolation points minus one to calculate distance between each
			//interpolation point and also total distance
			for (size_t i = 0; i < numberOfInterpolationPoints - 1; i++)
			{
				// calculate difference in x- and y- direction
				double differenceX = interpolationPoints[0][i + 1] - interpolationPoints[0][i];
				double differenceY = interpolationPoints[1][i + 1] - interpolationPoints[1][i];

				// calculate distance between each interpolation point by applying pythagorean theorem
				// and then taking the square root
				distance[i] = sqrt(sqrt(pow(differenceX, 2) + pow(differenceY, 2)));
				
				// sum up total distance in each iteration
				totalDistance += distance[i];

			}
	   
			//loop over all interpolation points - 1 to calculate parameter positions based on centripetal technique
			for (size_t i = 1; i < numberOfInterpolationPoints; i++)
			{
				t_bar[i] = t_bar[i - 1] + distance[i - 1] / totalDistance;
			}

			//return vector of parameter positions
			return { t_bar };
		}

		std::vector<double> knotVectorUsingAveraging(const std::vector<double>& parameterPositions,
			size_t polynomialDegree)
		{
			// polynomial degree, provided by user
			double p = polynomialDegree;

			// size of t_bar is equal to n
			size_t n = parameterPositions.size();

			// calculate size of knot vector based on n and p
			size_t m = n + polynomialDegree + 1;

			// check to make sure provided polynomial degree isn't too high relative to
			// number of parameter positions (and also # of interpolation points)
			if ((p + 1) * 2 > m)
			{
				throw std::runtime_error("Error. Please enter a lower polynomial degree");
			}

			// declare a vector of doubles to store knot Vector of size m
			std::vector<double> knotVector(m);

			// set left side of knotVector = 0. Number of zero entries depends on p
			for (size_t i = 0; i < polynomialDegree + 1; i++)
			{
				knotVector[i] = 0.0;
			}

			// calculate middle entries of knotVector using averaging technique
			// index i represents each entry of the knot vector
			for (size_t i = polynomialDegree + 1; i < n; i++)
			{
				// index j represents number of entries that must be averaged, based on p
				for (size_t j = 1; j < polynomialDegree + 1; j++)
				{
					knotVector[i] += (1 / p) * parameterPositions[i - p - 1 + j];
				}
			}

			// set right side of knotVector = 1. Number of zero entries depends on p
			for (size_t i = n; i < m; i++)
			{
				knotVector[i] = 1.0;
			}

			// return knotVector
			return { knotVector };


		}

	} // namespace splinekernel
} // namespace cie