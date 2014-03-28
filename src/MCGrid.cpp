#include "MCGrid.h"
#include "MC_LookupTable.h"

#include <math.h>

#include <iostream>

//------------------------------------------------------------------------------
// Constructors / Destructor
//------------------------------------------------------------------------------
MCGrid::MCGrid()
: _nbVertices(0),
  _resX(0),
  _resY(0),
  _resZ(0),
  _cubeSize(0.0)
{}

MCGrid::MCGrid(double cubeSize, const Vector3DF& volumeMin, const Vector3DF& volumeMax)
: _nbVertices(0),
  _resX(0),
  _resY(0),
  _resZ(0),
  _cubeSize(0.0)
{
	initGrid(cubeSize, volumeMin, volumeMax);
}

MCGrid::~MCGrid()
{}

//------------------------------------------------------------------------------
// Public functions
//------------------------------------------------------------------------------
void MCGrid::initGrid(double cubeSize, const Vector3DF& volumeMin, const Vector3DF& volumeMax)
{

	// Compute number of cell on each axis
	_resX = static_cast<int>(ceil((volumeMax.x-volumeMin.x)/cubeSize));
	_resY = static_cast<int>(ceil((volumeMax.y-volumeMin.y)/cubeSize));
	_resZ = static_cast<int>(ceil((volumeMax.z-volumeMin.z)/cubeSize));

	//std::cout << "_resX = " << _resX << std::endl;
	//std::cout << "_resY = " << _resX << std::endl;
	//std::cout << "_resZ = " << _resX << std::endl;

	// Compute grid dimensions
	_cubeSize = cubeSize;
	_volMin = volumeMin;

	//std::cout << "_cubeSize = " << _cubeSize << std::endl;
	//std::cout << "_volumeMin = " << _volMin.x <<  ", " << _volMin.y << ", " << _volMin.z << std::endl;

	_dimensions = Vector3DF(_resX, _resY, _resZ);
	_dimensions *= _cubeSize;

	//std::cout << "_dimensions = " << _dimensions.x <<  ", " << _dimensions.y << ", " << _dimensions.z << std::endl;

	// Init vertices data
	_nbVertices = _resX*_resY*_resZ;
	_vertices.resize(_nbVertices);
	for (int v=0; v<_nbVertices; ++v)
	{
		_vertices[v] = -1;	// By default, there's no vertex data
	}

	//std::cout << "_nbVertices = " << _nbVertices << std::endl;

}

void MCGrid::triangulate(Mesh&						mesh,
						 std::vector<Vector3DF>&	normals,
						 bool						computeNormals)
{
	std::vector<Mesh::Triangle>&	triangles = mesh.triangles();
	std::vector<Vector3DF>&			points = mesh.points();

	// Compute normals at every grid vertices
	if (computeNormals)
	{
		updateNormals();
	}

	// Iterate all cubes with data
	int nbVerticesData = _verticesData.size();
	//std::cout << "nbVerticesData = " << nbVerticesData << std::endl;
	for (int v=0; v<nbVerticesData; ++v)
	{
		MCVertex &vertex = _verticesData[v];

		/*std::cout << "(" << v << ") " << 
			"value= " << vertex.value << ", " <<
			"gridIndex=" << vertex.gridIndex << std::endl;*/

		// Make sure the vertex is associated with a cell
		unsigned int ix = getXIndex(vertex.gridIndex);
		unsigned int iy = getYIndex(vertex.gridIndex);
		unsigned int iz = getZIndex(vertex.gridIndex);
		//if (ix>=_resX) std::cout << "Depassement en x" << std::endl;
		//if (iy>=_resY) std::cout << "Depassement en y" << std::endl;
		//if (iz>=_resZ) std::cout << "Depassement en z" << std::endl;
		if ((ix < (_resX-1)) && (iy < (_resY-1)) && (iz < (_resZ-1)))
		{
			//std::cout << "is associated with a cell" << std::endl;
			// Get cube vertices indices
			// See Lorensen1987 for details on notation
			int v1 = _vertices[getGridIndex(ix,   iy,   iz)];
			int v2 = _vertices[getGridIndex(ix+1, iy,   iz)];
			int v3 = _vertices[getGridIndex(ix+1, iy+1, iz)];
			int v4 = _vertices[getGridIndex(ix  , iy+1, iz)];
			int v5 = _vertices[getGridIndex(ix  , iy  , iz+1)];
			int v6 = _vertices[getGridIndex(ix+1, iy  , iz+1)];
			int v7 = _vertices[getGridIndex(ix+1, iy+1, iz+1)];
			int v8 = _vertices[getGridIndex(ix  , iy+1, iz+1)];

            /*if (v < 10000)
            {
                std::clog << "v1 = " << v1 << std::endl;
                std::clog << "v2 = " << v2 << std::endl;
                std::clog << "v3 = " << v3 << std::endl;
                std::clog << "v4 = " << v4 << std::endl;
                std::clog << "v5 = " << v5 << std::endl;
                std::clog << "v6 = " << v6 << std::endl;
                std::clog << "v7 = " << v7 << std::endl;
                std::clog << "v8 = " << v8 << std::endl;
            }*/
/*
			std::cout << "res = " << _resX << ", " << _resY << ", " << _resZ << std::endl;
			std::cout << "ix, iy, iz= " << ix << ", " << iy << ", " << iz << std::endl;

			std::cout << "getGridIndex = " << getGridIndex(ix,   iy,   iz) << std::endl;
			std::cout << "getGridIndex = " << getGridIndex(ix+1,   iy,   iz) << std::endl;
			std::cout << "getGridIndex = " << getGridIndex(ix+1,   iy+1,   iz) << std::endl;
			std::cout << "getGridIndex = " << getGridIndex(ix,   iy+1,   iz) << std::endl;
			std::cout << "getGridIndex = " << getGridIndex(ix,   iy,   iz+1) << std::endl;
			std::cout << "getGridIndex = " << getGridIndex(ix+1,   iy,   iz+1) << std::endl;
			std::cout << "getGridIndex = " << getGridIndex(ix+1,   iy+1,   iz+1) << std::endl;
			std::cout << "getGridIndex = " << getGridIndex(ix,   iy+1,   iz+1) << std::endl;

			std::cout << "v1= " << v1 << std::endl;
			std::cout << "v2= " << v2 << std::endl;
			std::cout << "v3= " << v3 << std::endl;
			std::cout << "v4= " << v4 << std::endl;
			std::cout << "v5= " << v5 << std::endl;
			std::cout << "v6= " << v6 << std::endl;
			std::cout << "v7= " << v7 << std::endl;
			std::cout << "v8= " << v8 << std::endl;
*/
			// Make sure all vertices has data
			if (v1>=0 && v2>=0 && v3>=0 && v4>=0 && v5>=0 && v6>=0 && v7>=0 && v8>=0)
			{
				// Get cell vertices data
				MCVertex &vertex1 = _verticesData[v1];
				MCVertex &vertex2 = _verticesData[v2];
				MCVertex &vertex3 = _verticesData[v3];
				MCVertex &vertex4 = _verticesData[v4];
				MCVertex &vertex5 = _verticesData[v5];
				MCVertex &vertex6 = _verticesData[v6];
				MCVertex &vertex7 = _verticesData[v7];
				MCVertex &vertex8 = _verticesData[v8];

				// Build index in MC lookup table
				unsigned int cubeIndex = 0;
				if (vertex1.value<0.0) cubeIndex |= 1;
				if (vertex2.value<0.0) cubeIndex |= 2;
				if (vertex3.value<0.0) cubeIndex |= 4;
				if (vertex4.value<0.0) cubeIndex |= 8;
				if (vertex5.value<0.0) cubeIndex |= 16;
				if (vertex6.value<0.0) cubeIndex |= 32;
				if (vertex7.value<0.0) cubeIndex |= 64;
				if (vertex8.value<0.0) cubeIndex |= 128;

                /*if (v < 10000)
                {
                    std::cout << "vertex1.value = " << vertex1.value << std::endl;
                    std::cout << "vertex2.value = " << vertex2.value << std::endl;
                    std::cout << "vertex3.value = " << vertex3.value << std::endl;
                    std::cout << "vertex4.value = " << vertex4.value << std::endl;
                    std::cout << "vertex5.value = " << vertex5.value << std::endl;
                    std::cout << "vertex6.value = " << vertex6.value << std::endl;
                    std::cout << "vertex7.value = " << vertex7.value << std::endl;
                    std::cout << "vertex8.value = " << vertex8.value << std::endl;
                }*/

                //std::cout << "cubeIndex: " << cubeIndex << std::endl;

				// Make sure it is a surface cell
				if ((cubeIndex != 0) && (cubeIndex != 255))
				{
					int i=0;
					while (MC_LookupTable::trianglesList[cubeIndex][i] != -1)
					{
						// Compute/get triangle points on edges
						int p1, p2, p3;
						if (computeNormals)
						{
							p1 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
											  MC_LookupTable::trianglesList[cubeIndex][i+0]+1,
											  points,
											  normals);
							p2 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
											  MC_LookupTable::trianglesList[cubeIndex][i+1]+1,
											  points,
											  normals);
							p3 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
											  MC_LookupTable::trianglesList[cubeIndex][i+2]+1,
											  points,
											  normals);
						}
						else
						{

							p1 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
											  MC_LookupTable::trianglesList[cubeIndex][i+0]+1,
											  points);
							p2 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
											  MC_LookupTable::trianglesList[cubeIndex][i+1]+1,
											  points);
							p3 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
											  MC_LookupTable::trianglesList[cubeIndex][i+2]+1,
											  points);

						}

						// Create triangle
						triangles.push_back(Mesh::Triangle(p1, p2, p3));

						i += 3;
					}
				}
			}
		}
	}

    std::cout << "nbPoints: " << points.size() << std::endl;
}

//------------------------------------------------------------------------------
// Private functions
//------------------------------------------------------------------------------
void MCGrid::updateNormals()
{
	// Compute normals for every vertices with data
	for (long v=0; v<getNbVertices(); ++v)
	{
		// Only process vertices with data
		if (_vertices[v] == -1)
		{
			continue;
		}

		unsigned int ix = getXIndex(v);
		unsigned int iy = getYIndex(v);
		unsigned int iz = getZIndex(v);

		double value = getScalarValue(ix, iy, iz);

		Vector3DF& normal = getNormal(ix, iy, iz);

		// Compute df/dx
		// NOTE: We use a central difference to evaluate the gradient.
		//       Furthermore, we do not divide by dx, since we will normalize
		//       the normal afterward (this only works because dx==dy==dz).
		int previousDataIndex = (ix>0) ? _vertices[getGridIndex(ix-1, iy, iz)] : -1;
		int nextDataIndex = (ix<_resX-1) ? _vertices[getGridIndex(ix+1, iy, iz)] : -1;
		if (previousDataIndex == -1)
		{
			if (nextDataIndex == -1)
			{
				// If previous and next cells have no value
				// we consider that there is no variation
				normal.x = 0.0;
			}
			else
			{
				// Extrapolate fromt ix and ix+1
				normal.x = 2.0 * (_verticesData[nextDataIndex].value - value);
			}
		}
		else if (nextDataIndex == -1)
		{
			// Extrapolate from ix-1 and ix
			normal.x = 2.0 * (value - _verticesData[previousDataIndex].value);
		}
		else
		{
			normal.x = _verticesData[nextDataIndex].value -
				_verticesData[previousDataIndex].value;
		}

		// Compute df/dy
		previousDataIndex = (iy>0) ? _vertices[getGridIndex(ix, iy-1, iz)] : -1;
		nextDataIndex = (iy<_resY-1) ? _vertices[getGridIndex(ix, iy+1, iz)] : -1;
		if (previousDataIndex == -1)
		{
			if (nextDataIndex == -1)
			{
				// If previous and next cells have no value
				// we consider that there is no variation
				normal.y = 0.0;
			}
			else
			{
				// Extrapolate fromt iy and iy+1
				normal.y = 2.0 * (_verticesData[nextDataIndex].value - value);
			}
		}
		else if (nextDataIndex == -1)
		{
			// Extrapolate from iy-1 and iy
			normal.y = 2.0 * (value - _verticesData[previousDataIndex].value);
		}
		else
		{
			normal.y = _verticesData[nextDataIndex].value -
				_verticesData[previousDataIndex].value;
		}

		// Compute df/dz
		previousDataIndex = (iz>0) ? _vertices[getGridIndex(ix, iy, iz-1)] : -1;
		nextDataIndex = (iz<_resZ-1) ? _vertices[getGridIndex(ix, iy, iz+1)] : -1;
		if (previousDataIndex == -1)
		{
			if (nextDataIndex == -1)
			{
				// If previous and next cells have no value
				// we consider that there is no variation
				normal.z = 0.0;
			}
			else
			{
				// Extrapolate fromt iy and iy+1
				normal.z = 2.0 * (_verticesData[nextDataIndex].value - value);
			}
		}
		else if (nextDataIndex == -1)
		{
			// Extrapolate from iy-1 and iy
			normal.z = 2.0 * (value - _verticesData[previousDataIndex].value);
		}
		else
		{
			normal.z = _verticesData[nextDataIndex].value -
				_verticesData[previousDataIndex].value;
		}

		// Normalize the normal!
		double length = sqrt(normal.x*normal.x + normal.y*normal.y + normal.z*normal.z);
		if (length != 0.0)
		{
			normal /= length;
		}
	}
}

//HMC_LookupTable::trianglesList[index][noSommet];	// noSommet, jusqu'a ce que val = -1
