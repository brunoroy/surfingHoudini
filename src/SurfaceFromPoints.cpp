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
												std::vector<Vector3DF>&			normals,
							 					double							resolution,
							 					double							influenceRadius,
												double							particleRadius,
												bool							computeNormals)
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
	computeIsoValuesScatter(mcgrid, points, resolution, influenceRadius, particleRadius);

	// Triangulate surface!!!!
	PerfTimerBasicLinux triangulateTimer;
	triangulateTimer.start();
	mcgrid.triangulate(mesh, normals, computeNormals);
	triangulateTimer.pause();

	// Print timer
	funcTimer.pause();
	std::cout << ">> Timer_triangulate: " << triangulateTimer.getElapsedTime() << std::endl;
	std::cout << "Timer_createSurfaceFromPoints: " << funcTimer.getElapsedTime() << std::endl;
}

void SurfaceFromPoints::createSurfaceFromPointsFaster(const std::vector<Vector3DF>&	points,
							 						  Mesh&							mesh,
													  std::vector<Vector3DF>&		normals,
							 						  double						resolution,
							 						  double						influenceRadius,
													  double						particleRadius,
									   				  const std::vector<double>&	pointMasses,
									   				  const std::vector<double>&	pointDensities,
									   				  double						sphRadius,
													  bool							useColorField,
													  bool							computeNormals)
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
	double borderSize = particleRadius;
	borderSize += resolution;
	Vector3DF radiusVolume(borderSize, borderSize, borderSize);
	volMin -= radiusVolume;
	volMax += radiusVolume;

	// Init MC Grid
	MCGrid mcgrid(resolution, volMin, volMax);

	// Scalar field computations
	std::vector<unsigned int> surfacePoints;
	std::vector<unsigned int> surfaceVertices;

	SpatialGrid<PointInfo> spatialGrid(influenceRadius, volMin, volMax);

	buildSpatialGrid(spatialGrid, points);

	if (useColorField)
	{
		findSurfacePointsCF(points, spatialGrid, pointMasses, pointDensities, sphRadius, surfacePoints);
	}
	else
	{
		findSurfacePoints(points, spatialGrid, surfacePoints);
	}

	findSurfaceVertices(points, surfacePoints, mcgrid, particleRadius*2.1, surfaceVertices);
	computeIsoValues(mcgrid, spatialGrid, surfaceVertices, resolution, influenceRadius, particleRadius);

	// Triangulate surface!!!!
	PerfTimerBasicLinux triangulateTimer;
	triangulateTimer.start();
	mcgrid.triangulate(mesh, normals, computeNormals);
	triangulateTimer.pause();

	// Print timer
	funcTimer.pause();
	std::cout << ">> Timer_triangulate: " << triangulateTimer.getElapsedTime() << std::endl;
	std::cout << "Timer_createSurfaceFromPointsFaster: " << funcTimer.getElapsedTime() << std::endl;
}

//------------------------------------------------------------------------------
// Private functions
//------------------------------------------------------------------------------
void SurfaceFromPoints::computeIsoValuesScatter(MCGrid&						mcgrid,
										 		const std::vector<Vector3DF>&	points,
										 		double							resolution,
										 		double							influenceRadius,
										 		double							particleRadius)
{
	double influenceRadius2 = influenceRadius*influenceRadius;
	double influenceRadius6 = pow(influenceRadius, 6);

	// Init cells properties used to compute iso value
	std::vector<double> 			sumWj;			// SUM(Wj)
	std::vector<Eigen::Matrix3d>	sumRjGradWjT;	// SUM(rj*Gradient(Wj)')
	std::vector<Vector3DF>			sumGradWj;		// SUM(Gradient(Wj))
	std::vector<Vector3DF> 			sumRjWj;		// SUM(rj*Wj)

	int nbGridVertices = mcgrid.getNbVertices();
	sumWj.resize(nbGridVertices, 0.0);
	sumRjGradWjT.resize(nbGridVertices, Eigen::Matrix3d::Zero());
	sumGradWj.resize(nbGridVertices, Vector3DF(0.0,0.0,0.0));
	sumRjWj.resize(nbGridVertices, Vector3DF(0.0,0.0,0.0));

    //std::cout << "influenceRadius2: " << influenceRadius2 << std::endl;

	// Traverse points and add their contribution to nearby cells
	int nbPoints = points.size();
    //std::cout << "nbPoints: " << nbPoints << std::endl;
	for (int p=0; p<nbPoints; ++p)
	{
		// Get Nearby cells
		unsigned int minX, maxX, minY, maxY, minZ, maxZ;
		mcgrid.getCellsInRadius(points[p], influenceRadius, minX, maxX, minY, maxY, minZ, maxZ);

        int c = 0;

		// Process nearby cells
		Vector3DF vertexPos;
		for (int iz=minZ; iz<=maxZ; ++iz)
		{
			for (int iy=minY; iy<=maxY; ++iy)
			{
				for (int ix=minX; ix<=maxX; ++ix)
				{
					unsigned int cellIndex = mcgrid.getGridIndex(ix, iy, iz);

                    /*if (p == 0)
                        std::cout << "p = [" << points[p].x  << "," << points[p].y << "," << points[p].z << "]" << std::endl;

                    if (p == 0)
                        std::cout << "cellIndex: " << cellIndex << std::endl;*/

					mcgrid.getVertexPosition(ix, iy, iz, vertexPos);

                    /*if (p == 0)
                        std::cout << "vertexPos: [" << vertexPos.x << "," << vertexPos.y << "," << vertexPos.z << "]" << std::endl;*/

					// Is cell inside influence radius?
					Vector3DF delta(vertexPos);
					delta -= points[p];

                    //if (p == 0)
                    //    std::cout << "delta: [" << delta.x << "," << delta.y << "," << delta.z << "]" << std::endl;

					double dist2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;

                    /*if (p == 0)
                        std::clog << "dist2: " << dist2 << std::endl;*/

					if (dist2 < influenceRadius2)
					{
						// Compute kernel and it's gradient
						double dist = sqrt(dist2);
						double Wj = pow((1.0 - pow(dist/influenceRadius,2)), 3);
                        /*if (p == 0)
                            std::cout << "Wj = " << Wj << std::endl;*/
						Vector3DF gradWj(delta);
						gradWj *= -6.0*pow(influenceRadius2-dist2, 2) / influenceRadius6;

						// Update summation terms of cell
						sumWj[cellIndex] += Wj;

						sumRjGradWjT[cellIndex](0,0) += points[p].x*gradWj.x;
						sumRjGradWjT[cellIndex](1,0) += points[p].x*gradWj.y;
						sumRjGradWjT[cellIndex](2,0) += points[p].x*gradWj.z;
						sumRjGradWjT[cellIndex](0,1) += points[p].y*gradWj.x;
						sumRjGradWjT[cellIndex](1,1) += points[p].y*gradWj.y;
						sumRjGradWjT[cellIndex](2,1) += points[p].y*gradWj.z;
						sumRjGradWjT[cellIndex](0,2) += points[p].z*gradWj.x;
						sumRjGradWjT[cellIndex](1,2) += points[p].z*gradWj.y;
						sumRjGradWjT[cellIndex](2,2) += points[p].z*gradWj.z;

                        /*if (p == 0)
                        {
                            std::cout << "sumRjGradWjT" << std::endl;
                            std::cout << std::scientific;
                            std::cout << "[" << sumRjGradWjT[cellIndex](0,0) << " " << sumRjGradWjT[cellIndex](0,1) << " " << sumRjGradWjT[cellIndex](0,2) << "]" << std::endl;
                            std::cout << "[" << sumRjGradWjT[cellIndex](1,0) << " " << sumRjGradWjT[cellIndex](1,1) << " " << sumRjGradWjT[cellIndex](1,2) << "]" << std::endl;
                            std::cout << "[" << sumRjGradWjT[cellIndex](2,0) << " " << sumRjGradWjT[cellIndex](2,1) << " " << sumRjGradWjT[cellIndex](2,2) << "]" << std::endl;
                        }*/
						
						sumGradWj[cellIndex] += gradWj;
						
						sumRjWj[cellIndex].x += points[p].x*Wj;
						sumRjWj[cellIndex].y += points[p].y*Wj;
						sumRjWj[cellIndex].z += points[p].z*Wj;

                        c++;
					}

				}
			}
		}

        //if (p == 0)
            //std::cout << "test count: " << c << std::endl;
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
		if (sumWj[c] > 0.0)
		{
            /*if (c < 10000)
                std::cout << "cell = " << c << std::endl;*/

			mcgrid.getVertexPosition(ix, iy, iz, vertexPos);

			// Compute average position ( SUM(rj*Wj)/SUM(Wj) )
			Vector3DF averagePosition(sumRjWj[c]);
			averagePosition /= sumWj[c];

			// Compute the gradient of the average position
			// (1/SUM(Wj)) * SUM(rj*gradWj') - (1/SUM(Wj)^2) * SUM(gradWj) * SUM(rj*Wj)'
			Eigen::Matrix3d sumGradWjSumRjWjT;	// SUM(gradWj) * SUM(rj*Wj)'
			sumGradWjSumRjWjT(0,0) = sumGradWj[c].x*sumRjWj[c].x;
			sumGradWjSumRjWjT(0,1) = sumGradWj[c].x*sumRjWj[c].y;
			sumGradWjSumRjWjT(0,2) = sumGradWj[c].x*sumRjWj[c].z;
			sumGradWjSumRjWjT(1,0) = sumGradWj[c].y*sumRjWj[c].x;
			sumGradWjSumRjWjT(1,1) = sumGradWj[c].y*sumRjWj[c].y;
			sumGradWjSumRjWjT(1,2) = sumGradWj[c].y*sumRjWj[c].z;
			sumGradWjSumRjWjT(2,0) = sumGradWj[c].z*sumRjWj[c].x;
			sumGradWjSumRjWjT(2,1) = sumGradWj[c].z*sumRjWj[c].y;
			sumGradWjSumRjWjT(2,2) = sumGradWj[c].z*sumRjWj[c].z;

            /*if (c == 4370)
            {
                std::cout << "cell = " << c << std::endl;
                std::cout << "[" << sumGradWjSumRjWjT(0,0) << " " << sumGradWjSumRjWjT(0,1) << " " << sumGradWjSumRjWjT(0,2) << "]" << std::endl;
                std::cout << "[" << sumGradWjSumRjWjT(1,0) << " " << sumGradWjSumRjWjT(1,1) << " " << sumGradWjSumRjWjT(1,2) << "]" << std::endl;
                std::cout << "[" << sumGradWjSumRjWjT(2,0) << " " << sumGradWjSumRjWjT(2,1) << " " << sumGradWjSumRjWjT(2,2) << "]" << std::endl;
            }*/

            /*Eigen::Matrix3d gradAvgPosition =
                ((1.0/sumWj[c]) * sumRjGradWjT[c]) - ((1.0/(sumWj[c]*sumWj[c])) * sumGradWjSumRjWjT);*/

            double apTerm1 = 1.0/sumWj[c];
            //glm::dmat3x3 apTerm2 = apTerm1 * sumRjGradWjT[c];
            Eigen::Matrix3d apTerm2;
            for (int i = 0; i <= 2; ++i)
                for (int j = 0; j <= 2; ++j)
                    apTerm2(i,j) = apTerm1 * sumRjGradWjT[c](i,j);

            double apTerm3 = 1.0/(sumWj[c]*sumWj[c]);
            Eigen::Matrix3d gradAvgPosition;// = apTerm2 - (apTerm3 * sumGradWjSumRjWjT);
            for (int i = 0; i <= 2; ++i)
                for (int j = 0; j <= 2; ++j)
                {
                    double apTerm2IJ = apTerm2(i,j);
                    double sumIJ = sumGradWjSumRjWjT(i,j);
                    gradAvgPosition(i,j) = apTerm2IJ - (apTerm3 * sumIJ);
                    /*if (c == 4370)
                    {
                        std::clog << "apTerm2IJ: " << apTerm2IJ << std::endl;
                        std::clog << "sumIJ: " << sumIJ << std::endl;
                        std::clog << "arrayGradAvgPosition[" << i << "][" << j << "] = " << gradAvgPosition(i,j) << std::endl;
                    }*/
                }

            /*if (c == 4370)
            {
                std::cout << "apTerm1 = " << apTerm1 << std::endl;
                std::cout << "apTerm2" << std::endl;
                std::cout << std::scientific;
                std::cout << "[" << apTerm2(0,0) << " " << apTerm2(0,1) << " " << apTerm2(0,2) << "]" << std::endl;
                std::cout << "[" << apTerm2(1,0) << " " << apTerm2(1,1) << " " << apTerm2(1,2) << "]" << std::endl;
                std::cout << "[" << apTerm2(2,0) << " " << apTerm2(2,1) << " " << apTerm2(2,2) << "]" << std::endl;
                std::cout << "apTerm3 = " << apTerm1 << std::endl;

                std::cout << "sumWj[c] = " << sumWj[c] << std::endl;

                std::cout << "sumRjGradWjT" << std::endl;
                std::cout << std::scientific;
                std::cout << "[" << sumRjGradWjT[c](0,0) << " " << sumRjGradWjT[c](0,1) << " " << sumRjGradWjT[c](0,2) << "]" << std::endl;
                std::cout << "[" << sumRjGradWjT[c](1,0) << " " << sumRjGradWjT[c](1,1) << " " << sumRjGradWjT[c](1,2) << "]" << std::endl;
                std::cout << "[" << sumRjGradWjT[c](2,0) << " " << sumRjGradWjT[c](2,1) << " " << sumRjGradWjT[c](2,2) << "]" << std::endl;

                std::cout << "sumGradWjSumRjWjT" << std::endl;
                std::cout << std::scientific;
                std::cout << "[" << sumGradWjSumRjWjT(0,0) << " " << sumGradWjSumRjWjT(0,1) << " " << sumGradWjSumRjWjT(0,2) << "]" << std::endl;
                std::cout << "[" << sumGradWjSumRjWjT(1,0) << " " << sumGradWjSumRjWjT(1,1) << " " << sumGradWjSumRjWjT(1,2) << "]" << std::endl;
                std::cout << "[" << sumGradWjSumRjWjT(2,0) << " " << sumGradWjSumRjWjT(2,1) << " " << sumGradWjSumRjWjT(2,2) << "]" << std::endl;


                std::cout << "gradAvgPosition" << std::endl;
                std::cout << std::scientific;
                std::cout << "[" << gradAvgPosition(0,0) << " " << gradAvgPosition(0,1) << " " << gradAvgPosition(0,2) << "]" << std::endl;
                std::cout << "[" << gradAvgPosition(1,0) << " " << gradAvgPosition(1,1) << " " << gradAvgPosition(1,2) << "]" << std::endl;
                std::cout << "[" << gradAvgPosition(2,0) << " " << gradAvgPosition(2,1) << " " << gradAvgPosition(2,2) << "]" << std::endl;
            }*/

            /*if (c == 4370)
            {
                std::cout << "cell = " << c << std::endl;
                std::cout << "[" << gradAvgPosition(0,0) << " " << gradAvgPosition(0,1) << " " << gradAvgPosition(0,2) << "]" << std::endl;
                std::cout << "[" << gradAvgPosition(1,0) << " " << gradAvgPosition(1,1) << " " << gradAvgPosition(1,2) << "]" << std::endl;
                std::cout << "[" << gradAvgPosition(2,0) << " " << gradAvgPosition(2,1) << " " << gradAvgPosition(2,2) << "]" << std::endl;
            }*/

			// Find maximum eigenvalue of the gradient using the
			// Power method 
			double x[3] = { 1.0, 1.0, 1.0 };
			double newX[3];
			double maxValue, oldMaxValue = std::numeric_limits<double>::max();
			double threshold = 0.00001;
			double error = std::numeric_limits<double>::max();
			for (int i=0; (error > threshold) && i<500; ++i)
			{
				newX[0] = gradAvgPosition(0,0)*x[0] + gradAvgPosition(0,1)*x[1] + gradAvgPosition(0,2)*x[2];
				newX[1] = gradAvgPosition(1,0)*x[0] + gradAvgPosition(1,1)*x[1] + gradAvgPosition(1,2)*x[2];
				newX[2] = gradAvgPosition(2,0)*x[0] + gradAvgPosition(2,1)*x[1] + gradAvgPosition(2,2)*x[2];

				double absNewX0 = fabs(newX[0]);
				double absNewX1 = fabs(newX[1]);
				double absNewX2 = fabs(newX[2]);
				if ( (absNewX0 >= absNewX1) && (absNewX0 >= absNewX2) )
				{
					maxValue = newX[0];
				}
				else if (absNewX1 >= absNewX2)
				{
					maxValue = newX[1];
				}
				else
				{
					maxValue = newX[2];
				}

				if (maxValue==0.0)
				{
					break;
				}

				x[0] = newX[0] / maxValue;
				x[1] = newX[1] / maxValue;
				x[2] = newX[2] / maxValue;

				if (i>0)
				{
					error = fabs(maxValue-oldMaxValue);
					oldMaxValue = maxValue;
				}
				else
				{
					oldMaxValue = maxValue;
				}

                /*if (c == 0 && i == 499)
                {
                    std::cout << "x = [" << x[0] << "," << x[1] << "," << x[2] << "]" << std::endl;
                    std::cout << "i = " << i << std::endl;
                    std::cout << "error = " << error << std::endl;
                    std::cout << "threshold = " << threshold << std::endl;
                }*/

				// TODO: We could check if maxValue moves away from range [tlow, thigh] and
				// terminate earlier the algorithm! (Smarter, faster!)
			}

			double EVmax = fabs(maxValue);

			// Compute Radius correction based on EVmax
			double f = 1.0;
			const double tHigh = 2.0;
			const double tLow = 0.4;
			if (EVmax > tLow)
			{
				if (EVmax > tHigh) EVmax = tHigh;
				double gamma = (tHigh-EVmax) / (tHigh-tLow);
				f = gamma*gamma*gamma - 3.0*gamma*gamma + 3.0*gamma;
			}

			// Compute isoValue!!! (Finally...)
			Vector3DF deltaToAverage(vertexPos);
			deltaToAverage -= averagePosition;

			isoValue = sqrt(deltaToAverage.x*deltaToAverage.x + 
								   deltaToAverage.y*deltaToAverage.y +
								   deltaToAverage.z*deltaToAverage.z);
			isoValue -= particleRadius*f;
		}

		// Set value
		mcgrid.setScalarValue(ix, iy, iz, isoValue);
	}
}

void SurfaceFromPoints::computeIsoValues(MCGrid&					mcgrid,
										 SpatialGrid<PointInfo>&	spatialGrid,
						  				 std::vector<unsigned int>&	surfaceVertices,
										 double						resolution,
										 double						influenceRadius,
										 double						particleRadius)
{
	PerfTimerBasicLinux funcTimer;
	funcTimer.start();

	PerfTimerBasicLinux getElementsTimer;
	PerfTimerBasicLinux	contributionsTimer;
	PerfTimerBasicLinux	valuesTimer;

	double influenceRadius2 = influenceRadius*influenceRadius;
	double influenceRadius6 = pow(influenceRadius, 6);
	std::vector<PointInfo*>	nearbyPoints;

	// Init cells properties used to compute iso value
	std::vector<double> 			sumWj;			// SUM(Wj)
	std::vector<Eigen::Matrix3d>	sumRjGradWjT;	// SUM(rj*Gradient(Wj)')
	std::vector<Vector3DF>			sumGradWj;		// SUM(Gradient(Wj))
	std::vector<Vector3DF> 			sumRjWj;		// SUM(rj*Wj)

	int nbGridVertices = mcgrid.getNbVertices();
	sumWj.resize(nbGridVertices, 0.0);
	sumRjGradWjT.resize(nbGridVertices, Eigen::Matrix3d::Zero());
	sumGradWj.resize(nbGridVertices, Vector3DF(0.0,0.0,0.0));
	sumRjWj.resize(nbGridVertices, Vector3DF(0.0,0.0,0.0));

	// Compute cells isoValues
	Vector3DF vertexPos;
	for (int c=0; c<surfaceVertices.size(); ++c)
	{
		unsigned int cellIndex = surfaceVertices[c];
		unsigned int ix = mcgrid.getXIndex(cellIndex);
		unsigned int iy = mcgrid.getYIndex(cellIndex);
		unsigned int iz = mcgrid.getZIndex(cellIndex);
		mcgrid.getVertexPosition(ix, iy, iz, vertexPos);

		// Compute contribution of nearby particles
		getElementsTimer.start();
		spatialGrid.getElements(vertexPos, influenceRadius, nearbyPoints);
		std::vector<PointInfo*>::iterator it = nearbyPoints.begin();
		std::vector<PointInfo*>::iterator itEnd = nearbyPoints.end();
		getElementsTimer.pause();
		contributionsTimer.start();
		for ( ; it!=itEnd; ++it)
		{
			Vector3DF point = (*it)->pos;
			
			// Is cell inside influence radius?
			Vector3DF delta(vertexPos);
			delta -= point;

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

				sumRjGradWjT[cellIndex](0,0) += point.x*gradWj.x;
				sumRjGradWjT[cellIndex](1,0) += point.x*gradWj.y;
				sumRjGradWjT[cellIndex](2,0) += point.x*gradWj.z;
				sumRjGradWjT[cellIndex](0,1) += point.y*gradWj.x;
				sumRjGradWjT[cellIndex](1,1) += point.y*gradWj.y;
				sumRjGradWjT[cellIndex](2,1) += point.y*gradWj.z;
				sumRjGradWjT[cellIndex](0,2) += point.z*gradWj.x;
				sumRjGradWjT[cellIndex](1,2) += point.z*gradWj.y;
				sumRjGradWjT[cellIndex](2,2) += point.z*gradWj.z;
				
				sumGradWj[cellIndex] += gradWj;

				sumRjWj[cellIndex].x += point.x*Wj;
				sumRjWj[cellIndex].y += point.y*Wj;
				sumRjWj[cellIndex].z += point.z*Wj;
			}
		}
		contributionsTimer.pause();

		valuesTimer.start();

		// Compute isoValue based on particles contributions
		double isoValue = 1.0;
		if (sumWj[cellIndex] > 0.0)
		{
			mcgrid.getVertexPosition(ix, iy, iz, vertexPos);

			// Compute average position ( SUM(rj*Wj)/SUM(Wj) )
			Vector3DF averagePosition(sumRjWj[cellIndex]);
			averagePosition /= sumWj[cellIndex];

			// Compute the gradient of the average position
			// (1/SUM(Wj)) * SUM(rj*gradWj') - (1/SUM(Wj)^2) * SUM(gradWj) * SUM(rj*Wj)'
			Eigen::Matrix3d sumGradWjSumRjWjT;	// SUM(gradWj) * SUM(rj*Wj)'
			sumGradWjSumRjWjT(0,0) = sumGradWj[cellIndex].x*sumRjWj[cellIndex].x;
			sumGradWjSumRjWjT(0,1) = sumGradWj[cellIndex].x*sumRjWj[cellIndex].y;
			sumGradWjSumRjWjT(0,2) = sumGradWj[cellIndex].x*sumRjWj[cellIndex].z;
			sumGradWjSumRjWjT(1,0) = sumGradWj[cellIndex].y*sumRjWj[cellIndex].x;
			sumGradWjSumRjWjT(1,1) = sumGradWj[cellIndex].y*sumRjWj[cellIndex].y;
			sumGradWjSumRjWjT(1,2) = sumGradWj[cellIndex].y*sumRjWj[cellIndex].z;
			sumGradWjSumRjWjT(2,0) = sumGradWj[cellIndex].z*sumRjWj[cellIndex].x;
			sumGradWjSumRjWjT(2,1) = sumGradWj[cellIndex].z*sumRjWj[cellIndex].y;
			sumGradWjSumRjWjT(2,2) = sumGradWj[cellIndex].z*sumRjWj[cellIndex].z;

			Eigen::Matrix3d gradAvgPosition = 
				((1.0/sumWj[cellIndex]) * sumRjGradWjT[cellIndex])
				- ((1.0/(sumWj[cellIndex]*sumWj[cellIndex])) * sumGradWjSumRjWjT);

			// Find maximum eigenvalue of the gradient using the
			// Power method 
			double x[3] = { 1.0, 1.0, 1.0 };
			double newX[3];
			double maxValue, oldMaxValue = std::numeric_limits<double>::max();
			double threshold = 0.00001;
			double error = std::numeric_limits<double>::max();
			for (int i=0; (error > threshold) && i<500; ++i)
			{
				newX[0] = gradAvgPosition(0,0)*x[0] + gradAvgPosition(0,1)*x[1] + gradAvgPosition(0,2)*x[2];
				newX[1] = gradAvgPosition(1,0)*x[0] + gradAvgPosition(1,1)*x[1] + gradAvgPosition(1,2)*x[2];
				newX[2] = gradAvgPosition(2,0)*x[0] + gradAvgPosition(2,1)*x[1] + gradAvgPosition(2,2)*x[2];

				double absNewX0 = fabs(newX[0]);
				double absNewX1 = fabs(newX[1]);
				double absNewX2 = fabs(newX[2]);
				if ( (absNewX0 >= absNewX1) && (absNewX0 >= absNewX2) )
				{
					maxValue = newX[0];
				}
				else if (absNewX1 >= absNewX2)
				{
					maxValue = newX[1];
				}
				else
				{
					maxValue = newX[2];
				}

				if (maxValue==0.0)
				{
					break;
				}

				x[0] = newX[0] / maxValue;
				x[1] = newX[1] / maxValue;
				x[2] = newX[2] / maxValue;

				if (i>0)
				{
					error = fabs(maxValue-oldMaxValue);
					oldMaxValue = maxValue;
				}
				else
				{
					oldMaxValue = maxValue;
				}
				// TODO: We could check if maxValue moves away from range [tlow, thigh] and
				// terminate earlier the algorithm! (Smarter, faster!)
			}

			double EVmax = fabs(maxValue);

			// Compute Radius correction based on EVmax
			double f = 1.0;
			const double tHigh = 2.0;
			const double tLow = 0.4;
			if (EVmax > tLow)
			{
				if (EVmax > tHigh) EVmax = tHigh;
				double gamma = (tHigh-EVmax) / (tHigh-tLow);
				f = gamma*gamma*gamma - 3.0*gamma*gamma + 3.0*gamma;
			}

			// Compute isoValue!!! (Finally...)
			Vector3DF deltaToAverage(vertexPos);
			deltaToAverage -= averagePosition;

			isoValue = sqrt(deltaToAverage.x*deltaToAverage.x + 
								   deltaToAverage.y*deltaToAverage.y +
								   deltaToAverage.z*deltaToAverage.z);
			isoValue -= particleRadius*f;
		}

		// Set value
		mcgrid.setScalarValue(ix, iy, iz, isoValue);

		valuesTimer.pause();
	}

	// Print timer
	funcTimer.pause();
	std::cout << ">>>> Timer_getElements: " << getElementsTimer.getElapsedTime() << std::endl;
	std::cout << ">>>> Timer_contributions: " << contributionsTimer.getElapsedTime() << std::endl;
	std::cout << ">>>> Timer_values: " << valuesTimer.getElapsedTime() << std::endl;
	std::cout << ">> Timer_computeIsoValues: " << funcTimer.getElapsedTime() << std::endl;
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

void SurfaceFromPoints::findSurfacePointsCF(const std::vector<Vector3DF>&	points,
											SpatialGrid<PointInfo>&		spatialGrid,
											const std::vector<double>&	pointMasses,
											const std::vector<double>&	pointDensities,
											double						sphRadius,
											std::vector<unsigned int>&	surfacePoints)
{
	PerfTimerBasicLinux funcTimer;
	funcTimer.start();

	// This method finds surface points using the smoothed color field method (Muller2003).
	// Here, we make the assumption that the density and particles mass are always constant
	// and remove them from the equation.
	// This assumption will identify particles with non-uniform distribution as surface
	// particles in some cases, which could make the global approach a little less effective.
	
	std::vector<PointInfo*>	neighbors;

	double h = sphRadius;
	double h2 = h*h;
	double precomputedKernel = -45.0 / (3.14159265359 * pow(h,6));
	
	// Init data structures
	_colorGradientNorm.resize(points.size());
	surfacePoints.clear();

	// Process all points
	int nbPoints = points.size();
	for (int p=0; p<nbPoints; ++p)
	{
		// Compute gradient of color field using neighbors contributions
		Vector3DF gradCs(0,0,0);
		spatialGrid.getElements(points[p], sphRadius, neighbors);
		bool isOnSurface = false;
		if (neighbors.size() >= 25)
		{
			for (int n=0; n<neighbors.size(); ++n)
			{
				unsigned int j = neighbors[n]->id;

				Vector3DF delta(points[p]);
				delta -= neighbors[n]->pos;


				double dist2 = delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
				if (dist2 < h2 && dist2 > 0.0)
				{
					double dist = sqrt(dist2);

					delta *= (pointMasses[j]/pointDensities[j])*pow(h-dist,2) * precomputedKernel / dist;
					gradCs += delta;
				}
			}

			// Compute norm of gradient
			double norm = sqrt(gradCs.x*gradCs.x + gradCs.y*gradCs.y + gradCs.z*gradCs.z);
			_colorGradientNorm[p] = norm;

			// Test against threshold
			if (norm >= 4.0)
			{
				isOnSurface = true;
			}
		}
		else
		{
			isOnSurface = true;
		}

		if (isOnSurface)
		{
			surfacePoints.push_back(p);
		}
	}

	// Print timer
	funcTimer.pause();
	std::cout << ">> Timer_findSurfacePointsCF: " << funcTimer.getElapsedTime() << std::endl;
}

void SurfaceFromPoints::findSurfacePoints(const std::vector<Vector3DF>&	points,
										  SpatialGrid<PointInfo>&		spatialGrid,
										  std::vector<unsigned int>&	surfacePoints)
{
	PerfTimerBasicLinux funcTimer;
	funcTimer.start();

	std::vector<PointInfo*> elements;
	surfacePoints.clear();

	_surfacePoints.clear();
	_surfacePoints.resize(points.size(), 0.0);

	std::cout << "Res: " << spatialGrid.getResX() << ", " << spatialGrid.getResY() << ", " 
		<< spatialGrid.getResZ() << std::endl;
	std::cout << "Cellsize: " << spatialGrid.getCellSize() << std::endl;
	std::cout << "Volume start: " << spatialGrid.getVolumeStart().x << ", " << spatialGrid.getVolumeStart().y
		<< ", " << spatialGrid.getVolumeStart().z << std::endl;

	// Find non-empty cells (in spatial grid) with empty neighbors
	int resX = spatialGrid.getResX();
	int resY = spatialGrid.getResY();
	int resZ = spatialGrid.getResZ();
	for (int ix=0; ix<resX; ++ix)
	{
		for (int iy=0; iy<resY; ++iy)
		{
			for (int iz=0; iz<resZ; ++iz)
			{
				if (spatialGrid.isCellEmpty(ix, iy, iz))
					continue;

				// Check neighbors
				bool hasEmptyNeighbor = false;
				if (ix==0 || ix==resX-1 || iy==0 || iy==resY-1 || iz==0 || iz==resZ-1)
				{
					hasEmptyNeighbor = true;
				}
				else
				{
					for (int dx=-1; !hasEmptyNeighbor && dx<=1; ++dx)
					{
						for (int dy=-1; !hasEmptyNeighbor && dy<=1; ++dy)
						{
							for (int dz=-1; !hasEmptyNeighbor && dz<=1; ++dz)
							{
								if (dx==0 && dy==0 && dz==0)
									continue;

								if (spatialGrid.isCellEmpty(ix+dx, iy+dy, iz+dz))
								{
									hasEmptyNeighbor = true;
								}
							}
						}
					}
				}

				// Add cell's points to surface points list if it's a surface cell
				if (hasEmptyNeighbor)
				{
					spatialGrid.getElements(ix, iy, iz, elements);
					for (int e=0; e<elements.size(); ++e)
					{
						surfacePoints.push_back(elements[e]->id);

						_surfacePoints[elements[e]->id] = 1.0;
					}
				}
			}
		}
	}

	// Print timer
	funcTimer.pause();
	std::cout << ">> Timer_findSurfacePoints: " << funcTimer.getElapsedTime() << std::endl;
}

void SurfaceFromPoints::findSurfaceVertices(const std::vector<Vector3DF>&		points,
							 				const std::vector<unsigned int>&	surfacePoints,
							 				const MCGrid&						mcgrid,
											double								aabbLength,
							 				std::vector<unsigned int>&			surfaceVertices)
{
	PerfTimerBasicLinux funcTimer;
	funcTimer.start();

	_surfaceVertices.clear();

	// Init data structure
	surfaceVertices.clear();
	_isNearSurfaceParticle.resize(mcgrid.getNbVertices());
	for (int v=0; v<_isNearSurfaceParticle.size(); ++v)
	{
		_isNearSurfaceParticle[v] = false;
	}

	// Iterates surface points and find vertices nearby
	for (int i=0; i<surfacePoints.size(); ++i)
	{
		Vector3DF point = points[surfacePoints[i]];

		// Find vertices in aabb around surface particle
		unsigned int minX, maxX, minY, maxY, minZ, maxZ;
		mcgrid.getCellsInRadius(point, aabbLength, minX, maxX, minY, maxY, minZ, maxZ);

		// Mark them as surface vertices
		for (int iz=minZ; iz<=maxZ; ++iz)
		{
			for (int iy=minY; iy<=maxY; ++iy)
			{
				for (int ix=minX; ix<=maxX; ++ix)
				{
					unsigned int cellIndex = mcgrid.getGridIndex(ix, iy, iz);

					_isNearSurfaceParticle[cellIndex] = true;
				}
			}
		}
	}


	// Add marked particles to list
	for (int v=0; v<_isNearSurfaceParticle.size(); ++v)
	{
		if (_isNearSurfaceParticle[v])
		{
			surfaceVertices.push_back(v);

			unsigned int ix = mcgrid.getXIndex(v);
			unsigned int iy = mcgrid.getYIndex(v);
			unsigned int iz = mcgrid.getZIndex(v);

			Vector3DF pos;
			mcgrid.getVertexPosition(ix, iy, iz, pos);


			_surfaceVertices.push_back(pos);
		
		}
	}

	// Print timer
	funcTimer.pause();
	std::cout << ">> Timer_findSurfaceVertices: " << funcTimer.getElapsedTime() << std::endl;
}

