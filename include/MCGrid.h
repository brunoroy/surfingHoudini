#ifndef MCGRID_H
#define MCGRID_H

// NOTE: Some functions might be accessed quite frequently, therefore they have been inlined and
// 		 their definition can be found in MCGrid.hpp instead of MCGrid.cpp. This is because inlined
// 		 functions must be defined inside the header file (MCGrid.hpp is included at the end of this
// 		 header file)

#include <Math/vector.h>
using Math::Vector3DF;

#include "Mesh.h"

// Class MCGrid : Marching Cube Grid
class MCGrid
{
public:
	MCGrid();
	MCGrid(double cubeSize, const Vector3DF& volumeMin, const Vector3DF& volumeMax);
	~MCGrid();

	void initGrid(double cubeSize, const Vector3DF& volumeMin, const Vector3DF& volumeMax);

	inline int getNbVertices() const { return _resX*_resY*_resZ; }
	inline unsigned int getResX() const { return _resX; }
	inline unsigned int getResY() const { return _resY; }
	inline unsigned int getResZ() const { return _resZ; }
	inline double getCubeSize() const { return _cubeSize; }
	inline Vector3DF getVolumeStart() const { return _volMin; }
	inline Vector3DF getDimensions() const { return _dimensions; }

	void getCellsInRadius(const Vector3DF&	position,
						  double			radius,
						  unsigned int&		minX,
						  unsigned int&		maxX,
						  unsigned int&		minY,
						  unsigned int&		maxY,
						  unsigned int&		minZ,
						  unsigned int&		maxZ) const;	// See MCGrid.hpp

	void getVertexPosition(unsigned int ix,
						   unsigned int iy,
						   unsigned int iz,
						   Vector3DF& position) const;	// See MCGrid.hpp

	double getScalarValue(unsigned int ix, unsigned int iy, unsigned int iz) const;	// See MCGrid.hpp
	void setScalarValue(unsigned int ix, unsigned int iy, unsigned int iz, double value);		// See MCGrid.hpp

	Vector3DF& getNormal(unsigned int ix, unsigned int iy, unsigned int iz);	// See MCGrid.hpp
	const Vector3DF& getNormal(unsigned int ix, unsigned int iy, unsigned int iz) const;	// See MCGrid.hpp
	

    void triangulate(Mesh& mesh);

	unsigned int getGridIndex(unsigned int ix, unsigned int iy, unsigned int iz) const;	// See MCGrid.hpp
	unsigned int getXIndex(unsigned int gridIndex) const; // See MCGrid.hpp
	unsigned int getYIndex(unsigned int gridIndex) const; // See MCGrid.hpp
	unsigned int getZIndex(unsigned int gridIndex) const; // See MCGrid.hpp

private:
	struct MCVertex;


	unsigned int getEdgePoint(MCVertex& v1, MCVertex& v2, MCVertex& v3, MCVertex& v4,
							  MCVertex& v5, MCVertex& v6, MCVertex& v7, MCVertex& v8,
							  int edgeNo,
							  std::vector<Vector3DF>& points) const;

	unsigned int getEdgePoint(MCVertex& v1, MCVertex& v2, MCVertex& v3, MCVertex& v4,
							  MCVertex& v5, MCVertex& v6, MCVertex& v7, MCVertex& v8,
							  int edgeNo,
							  std::vector<Vector3DF>& points,
							  std::vector<Vector3DF>& normals) const;

	void updateNormals(void);

private:
	unsigned int 			_nbVertices;
	std::vector<int>		_vertices;
	std::vector<MCVertex>	_verticesData;

	unsigned int 	_resX;
	unsigned int 	_resY;
	unsigned int 	_resZ;
	double			_cubeSize;
	Vector3DF		_volMin;
	Vector3DF		_dimensions;

};

// Inline function body
#include "MCGrid.hpp"

#endif	// MCGRID_H
