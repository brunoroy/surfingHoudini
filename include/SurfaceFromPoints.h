#ifndef SURFACEFROMPOINTS_H
#define SURFACEFROMPOINTS_H

#include "Mesh.h"

#include <Math/vector.h>
using Math::Vector3DF;

#include <vector>

class MCGrid;
template<class T> class SpatialGrid;

class SurfaceFromPoints
{
public:
	SurfaceFromPoints();
	~SurfaceFromPoints();

	void createSurfaceFromPoints(const std::vector<Vector3DF>&	points,
								 Mesh&							mesh,
                                 double 						resolution);

private:
	struct PointInfo;

	void computeIsoValuesScatter(MCGrid&						mcgrid,
						  		const	std::vector<Vector3DF>&	points,
                                double							resolution);

	void buildSpatialGrid(SpatialGrid<PointInfo>&		spatialGrid,
						  const std::vector<Vector3DF>&	points);


public:
	std::vector<Vector3DF>	_surfaceVertices;
	std::vector<double>		_surfacePoints;


private:
	std::vector<bool>	_isNearSurfaceParticle;
};


#endif	// SURFACEFROMPOINTS_H
