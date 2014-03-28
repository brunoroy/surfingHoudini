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
            isoValue -= resolution*f;
		}

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

