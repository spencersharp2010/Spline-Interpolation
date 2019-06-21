#include "basisfunctions.hpp"
#include "surface.hpp"

namespace cie
{
namespace splinekernel
{
namespace detail
{

std::vector<std::vector<double>> evaluateShapeFunctions( const std::vector<double>& knotVector,
                                                         size_t numberOfControlPoints ,
                                                         size_t numberOfSamplePoints )
{
    std::vector<std::vector<double>> shapeFunctionValues( numberOfSamplePoints );

    size_t polynomialDegree = knotVector.size( ) - numberOfControlPoints - 1;

    for( size_t iEvaluationCoordinate = 0; iEvaluationCoordinate < numberOfSamplePoints; ++iEvaluationCoordinate )
    {
        shapeFunctionValues[iEvaluationCoordinate].resize( numberOfControlPoints );

        double t = iEvaluationCoordinate / ( numberOfSamplePoints - 1.0 );

        for( size_t iShapeFunction = 0; iShapeFunction < numberOfControlPoints; ++iShapeFunction )
        {
            double value = evaluateBSplineBasis( t, iShapeFunction, polynomialDegree, knotVector );

            shapeFunctionValues[iEvaluationCoordinate][iShapeFunction] = value;
        }
    }

    return shapeFunctionValues;
}

double computeComponent( const std::vector<double>& Nr,
                         const std::vector<double>& Ns,
                         const linalg::Matrix& controlPointValues )
{
    size_t size1 = controlPointValues.size1( );
    size_t size2 = controlPointValues.size2( );

    double value = 0.0;

    // Compute sum over basis functions (tensor product) times control point value
    for( size_t iCP = 0; iCP < size1; ++iCP )
    {
        for( size_t jCP = 0; jCP < size2; ++jCP )
        {
            value += Nr[iCP] * Ns[jCP] * controlPointValues( iCP, jCP );
        }
    }

    return value;
}

} // splinesurfacehelper

VectorOfMatrices evaluateSurface( const std::array<std::vector<double>, 2>& knotVectors,
                                  const VectorOfMatrices& controlPoints,
                                  std::array<size_t, 2> numberOfSamplePoints )
{
    using Shapes1D = std::vector<std::vector<double>>;

    // First evaluate shape functions separately in both coordinate directions
    Shapes1D shapesR = detail::evaluateShapeFunctions( knotVectors[0], controlPoints[0].size1( ), numberOfSamplePoints[0] );
    Shapes1D shapesS = detail::evaluateShapeFunctions( knotVectors[1], controlPoints[0].size2( ), numberOfSamplePoints[1] );

    VectorOfMatrices result( controlPoints.size( ) );

    // Loop over components, e.g. x, y and z, each being a 2D matrix of values
    for( size_t iComponent = 0; iComponent < controlPoints.size( ); ++iComponent )
    {
        result[iComponent] = linalg::Matrix( numberOfSamplePoints[0], numberOfSamplePoints[1], 0.0 );

        // Loop over all sample points in local coordinates r and s
        for( size_t iR = 0; iR < numberOfSamplePoints[0]; ++iR )
        {
            for( size_t iS = 0; iS < numberOfSamplePoints[1]; ++iS )
            {
                // Compute tensor product and multiply by control point values
                result[iComponent]( iR, iS ) = detail::computeComponent( shapesR[iR], shapesS[iS], controlPoints[iComponent] );
            }
        }
    }

    return result;
}

} // namespace splinekernel
} // namespace cie
