#include "interpolation.hpp"
#include "basisfunctions.hpp"
#include "linalg.hpp"

#include <cmath>
#include <exception>

namespace cie
{
namespace splinekernel
{

ControlPointsAndKnotVector interpolateWithBSplineCurve( const ControlPoints2D& interpolationPoints,
                                                        size_t polynomialDegree )
{
    // Throw exception if number of x-values is not equal number of y-values
    
    return { { }, { } };
}

std::vector<double> centripetalParameterPositions( const ControlPoints2D& interpolationPoints )
{
    return { };
}

std::vector<double> knotVectorUsingAveraging( const std::vector<double>& parameterPositions,
                                              size_t polynomialDegree )
{
    // Throw exception if polynomial degree is too high for given number of points

    return { };
}

} // namespace splinekernel
} // namespace cie
