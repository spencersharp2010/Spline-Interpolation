#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"

#include "basisfunctions.hpp"
#include "curve.hpp"
#include "surface.hpp"
#include "interpolation.hpp"

// This header defines how to convert between numpy array and linalg::Matrix
#include "matrixConversion.hpp"

PYBIND11_MODULE( pysplinekernel, m ) 
{
    m.doc( ) = "spline computation kernel"; // optional module docstring

    m.def( "evaluateBSplineBasis", &cie::splinekernel::evaluateBSplineBasis, "Evaluates single b-spline basis function." );
    m.def( "evaluate2DCurve", &cie::splinekernel::evaluate2DCurve, "Evaluates B-Spline curve by multiplying control points and basis functions." );
    m.def( "evaluate2DCurveDeBoor", &cie::splinekernel::evaluate2DCurveDeBoor, "Evaluates B-Spline using DeBoor" );
    m.def( "evaluateSurface", &cie::splinekernel::evaluateSurface, "Evaluates B-Spline surface" );
	m.def( "evaluateInterpolateWithBSplineCurve", &cie::splinekernel::interpolateWithBSplineCurve, "Interpolates with B-Spline Curve");

}
