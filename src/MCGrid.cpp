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

	// Compute grid dimensions
	_cubeSize = cubeSize;
	_volMin = volumeMin;

	_dimensions = Vector3DF(_resX, _resY, _resZ);
	_dimensions *= _cubeSize;

	// Init vertices data
    _nbVertices = _resX*_resY*_resZ;
    _vertices.resize(_nbVertices);
	for (int v=0; v<_nbVertices; ++v)
	{
		_vertices[v] = -1;	// By default, there's no vertex data
	}

	//std::cout << "_nbVertices = " << _nbVertices << std::endl;

}

void MCGrid::triangulate(Mesh& mesh)
{
	std::vector<Mesh::Triangle>&	triangles = mesh.triangles();
	std::vector<Vector3DF>&			points = mesh.points();

	// Iterate all cubes with data
	int nbVerticesData = _verticesData.size();
	//std::cout << "nbVerticesData = " << nbVerticesData << std::endl;
    for (int v=0; v<nbVerticesData; ++v)
	{
		MCVertex &vertex = _verticesData[v];

		// Make sure the vertex is associated with a cell
		unsigned int ix = getXIndex(vertex.gridIndex);
		unsigned int iy = getYIndex(vertex.gridIndex);
		unsigned int iz = getZIndex(vertex.gridIndex);
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

				// Make sure it is a surface cell
                if ((cubeIndex != 0) && (cubeIndex != 255))
				{
					int i=0;
                    while (MC_LookupTable::trianglesList[cubeIndex][i] != -1)
                    {
						// Compute/get triangle points on edges
						int p1, p2, p3;
                        p1 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
                                          MC_LookupTable::trianglesList[cubeIndex][i+0]+1,
                                          points);
                        p2 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
                                          MC_LookupTable::trianglesList[cubeIndex][i+1]+1,
                                          points);
                        p3 = getEdgePoint(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8,
                                          MC_LookupTable::trianglesList[cubeIndex][i+2]+1,
                                          points);

						// Create triangle
                        triangles.push_back(Mesh::Triangle(p1, p2, p3));

                        i += 3;
                    }
				}
            }
		}
	}
}

//HMC_LookupTable::trianglesList[index][noSommet];	// noSommet, jusqu'a ce que val = -1
