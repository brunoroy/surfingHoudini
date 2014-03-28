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
								 std::vector<Vector3DF>&		normals,
								 double 						resolution,
								 double							influenceRadius,
								 double							particleRadius,
								 bool							computeNormals);

	void createSurfaceFromPointsFaster(const std::vector<Vector3DF>&	points,
									   Mesh&							mesh,
									   std::vector<Vector3DF>&			normals,
									   double							resolution,
									   double							influenceRadius,
									   double							particleRadius,
									   const std::vector<double>&		pointMasses,
									   const std::vector<double>&		pointDensities,
									   double							sphRadius,
									   bool								useColorField,
									   bool								computeNormals);

private:
	struct PointInfo;

	void computeIsoValuesScatter(MCGrid&						mcgrid,
						  		const	std::vector<Vector3DF>&	points,
						  		double							resolution,
						  		double							influenceRadius,
						  		double							particleRadius);

	void computeIsoValues(MCGrid&						mcgrid,
						  SpatialGrid<PointInfo>&		spatialGrid,
						  std::vector<unsigned int>&	surfaceVertices,
						  double						resolution,
						  double						influenceRadius,
						  double						particleRadius);

	void buildSpatialGrid(SpatialGrid<PointInfo>&		spatialGrid,
						  const std::vector<Vector3DF>&	points);

	void findSurfacePointsCF(const std::vector<Vector3DF>&	points,
							 SpatialGrid<PointInfo>&		spatialGrid,
							 const std::vector<double>&		pointMasses,
							 const std::vector<double>&		pointDensities,
							 double							sphRadius,
							 std::vector<unsigned int>&		surfacePoints);

	void findSurfacePoints(const std::vector<Vector3DF>&	points,
						   SpatialGrid<PointInfo>&			spatialGrid,
						   std::vector<unsigned int>&		surfacePoints);

	void findSurfaceVertices(const std::vector<Vector3DF>&		points,
							 const std::vector<unsigned int>&	surfacePoints,
							 const MCGrid&						mcgrid,
							 double								aabbLength,
							 std::vector<unsigned int>&			surfaceVertices);


public:
	std::vector<double>		_colorGradientNorm;
	std::vector<Vector3DF>	_surfaceVertices;
	std::vector<double>		_surfacePoints;


private:
	std::vector<bool>	_isNearSurfaceParticle;
};


#endif	// SURFACEFROMPOINTS_H
