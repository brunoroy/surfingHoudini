#include "SurfaceFromPoints.h"
#include "MCGrid.h"
#include "SpatialGrid.h"

#include <Eigen/Dense>

#include "PerfTimerBasicLinux.h"

// Structure used to store points info in spatial grid
struct SurfaceFromPoints::PointInfo
{
	PointInfo(const Vector3DF& paramPos, unsigned int paramID) : pos(paramPos), id(paramID) {}

	Vector3DF		pos;
	unsigned int	id;
};

namespace
{
	inline
	double max(double a, double b)
	{
		return (a>b) ? a : b;
	}
};



//------------------------------------------------------------------------------
// Constructor / Destructor
//------------------------------------------------------------------------------
SurfaceFromPoints::SurfaceFromPoints()
{
}

SurfaceFromPoints::~SurfaceFromPoints()
{
}

//------------------------------------------------------------------------------
// Public functions
//------------------------------------------------------------------------------
void SurfaceFromPoints::createSurfaceFromPoints(const std::vector<Vector3DF>&	points,
							 					Mesh&							mesh,
                                                double							resolution)
{
	PerfTimerBasicLinux funcTimer;
	funcTimer.start();

	// Make sure we have data to process
	if (points.size() <=0)
	{
		return;
	}

	// Get boundary informations of volume containing points
	Vector3DF volMin = points[0];
	Vector3DF volMax = points[0];
	int nbPoints = points.size();
	for (int p=0; p<nbPoints; ++p)
	{
		Vector3DF point(points[p]);

		if (volMin.x>point.x)
			volMin.x = point.x;
		if (volMin.y>point.y)
			volMin.y = point.y;
		if (volMin.z>point.z)
			volMin.z = point.z;

		if (volMax.x<point.x)
			volMax.x = point.x;
		if (volMax.y<point.y)
			volMax.y = point.y;
		if (volMax.z<point.z)
			volMax.z = point.z;
	}

	// Expand volume by the influence radius
    double influenceRadius = resolution * 4.0;
    double borderSize = influenceRadius;
	borderSize += resolution;
    Vector3DF radiusVolume(borderSize, borderSize, borderSize);
	volMin -= radiusVolume;
    volMax += radiusVolume;

    std::clog << "volume: minimum[" << volMin.x << "," << volMin.y << "," << volMin.z << "]";
    std::clog << "; maximum[" << volMax.x << "," << volMax.y << "," << volMax.z << "]";

	// Init MC Grid
	MCGrid mcgrid(resolution, volMin, volMax);

	// Scalar field computations
    computeIsoValuesScatter(mcgrid, points, resolution);

	// Triangulate surface!!!!
	PerfTimerBasicLinux triangulateTimer;
	triangulateTimer.start();
    mcgrid.triangulate(mesh);
	triangulateTimer.pause();

	// Print timer
	funcTimer.pause();
	std::cout << ">> Timer_triangulate: " << triangulateTimer.getElapsedTime() << std::endl;
	std::cout << "Timer_createSurfaceFromPoints: " << funcTimer.getElapsedTime() << std::endl;
}

//------------------------------------------------------------------------------
// Private functions
//------------------------------------------------------------------------------
void SurfaceFromPoints::computeIsoValuesScatter(MCGrid&						mcgrid,
										 		const std::vector<Vector3DF>&	points,
                                                double							resolution)
{
    double influenceRadius = resolution * 4.0;
    double influenceRadius2 = influenceRadius*influenceRadius;
	double influenceRadius6 = pow(influenceRadius, 6);

	// Init cells properties used to compute iso value
	std::vector<double> 			sumWj;			// SUM(Wj)
    //std::vector<Eigen::Matrix3d>	sumRjGradWjT;	// SUM(rj*Gradient(Wj)')
    //std::vector<Vector3DF>			sumGradWj;		// SUM(Gradient(Wj))
    std::vector<Vector3DF> 			sumRjWj;		// SUM(rj*Wj)

    int nbGridVertices = mcgrid.getNbVertices();
    sumWj.resize(nbGridVertices, 0.0);
    //sumRjGradWjT.resize(nbGridVertices, Eigen::Matrix3d::Zero());
    //sumGradWj.resize(nbGridVertices, Vector3DF(0.0,0.0,0.0));
    sumRjWj.resize(nbGridVertices, Vector3DF(0.0,0.0,0.0));

	// Traverse points and add their contribution to nearby cells
	int nbPoints = points.size();
    std::cout << "nbPoints: " << nbPoints << std::endl;
    for (int p=0; p<nbPoints; ++p)
	{
		// Get Nearby cells
		unsigned int minX, maxX, minY, maxY, minZ, maxZ;
		mcgrid.getCellsInRadius(points[p], influenceRadius, minX, maxX, minY, maxY, minZ, maxZ);

        //int c = 0;

		// Process nearby cells
        Vector3DF vertexPos;
        for (int iz=minZ; iz<=maxZ; ++iz)
		{
			for (int iy=minY; iy<=maxY; ++iy)
			{
				for (int ix=minX; ix<=maxX; ++ix)
                {
					unsigned int cellIndex = mcgrid.getGridIndex(ix, iy, iz);

					mcgrid.getVertexPosition(ix, iy, iz, vertexPos);
					Vector3DF delta(vertexPos);
					delta -= points[p];

					double dist2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;

                    if (dist2 < influenceRadius2)
					{
						// Compute kernel and it's gradient
						double dist = sqrt(dist2);
						double Wj = pow((1.0 - pow(dist/influenceRadius,2)), 3);
						Vector3DF gradWj(delta);
						gradWj *= -6.0*pow(influenceRadius2-dist2, 2) / influenceRadius6;

						// Update summation terms of cell
						sumWj[cellIndex] += Wj;
						
						sumRjWj[cellIndex].x += points[p].x*Wj;
						sumRjWj[cellIndex].y += points[p].y*Wj;
                        sumRjWj[cellIndex].z += points[p].z*Wj;
                    }

				}
            }
        }
    }

	// Compute cells isoValues
	Vector3DF vertexPos;
    for (int c=0; c<nbGridVertices; ++c)
	{
        unsigned int ix = mcgrid.getXIndex(c);
        unsigned int iy = mcgrid.getYIndex(c);
        unsigned int iz = mcgrid.getZIndex(c);

        // Make sure there was contribution from at least one particle
        double isoValue = 1.0;
        /*if (sumWj[c] > 0.0)
        {*/
            mcgrid.getVertexPosition(ix, iy, iz, vertexPos);

            // Compute average position ( SUM(rj*Wj)/SUM(Wj) )
            Vector3DF averagePosition(sumRjWj[c]);
            averagePosition /= sumWj[c];

            // Compute isoValue!!! (Finally...)
            Vector3DF deltaToAverage(vertexPos);
            deltaToAverage -= averagePosition;

            isoValue = sqrt(deltaToAverage.x*deltaToAverage.x +
                                   deltaToAverage.y*deltaToAverage.y +
                                   deltaToAverage.z*deltaToAverage.z);
            isoValue -= resolution;//*f;
        //}

        // Set value
        mcgrid.setScalarValue(ix, iy, iz, isoValue);
	}
}

void SurfaceFromPoints::buildSpatialGrid(SpatialGrid<PointInfo>&		spatialGrid,
										 const std::vector<Vector3DF>&	points)
{
	PerfTimerBasicLinux funcTimer;
	funcTimer.start();

	// Iterate over all points and add them in the spatial grid
	int nbPoints = points.size();
	for (int p=0; p<nbPoints; ++p)
	{
		spatialGrid.insert(PointInfo(points[p], p), points[p]);
	}

	// Print timer
	funcTimer.pause();
	std::cout << ">> Timer_buildSpatialGrid: " << funcTimer.getElapsedTime() << std::endl;
}

