#ifndef SPATIALGRID_H
#define SPATIALGRID_H

// NOTE: Due to inlining and use of templates, all the function definitions can be found
// 		 in SpatialGrid.hpp which is "included" at the end of the file.

#include <Math/vector.h>
using Math::Vector3DF;

#include <list>
#include <vector>

// Class SpatialGrid : Grid for spatial partitionning
template<class T>
class SpatialGrid
{
public:
	SpatialGrid(double cellSize, const Vector3DF& volumeMin, const Vector3DF& volumeMax);
	~SpatialGrid();

	void initGrid(double cellSize, const Vector3DF& volumeMin, const Vector3DF& volumeMax);
	void clear();

	void insert(const T& element, const Vector3DF& position);
	void insert(const T& element, const Vector3DF& AABBmin, const Vector3DF& AABBmax);

	void getElements(const Vector3DF& position, double radius, std::vector<T*>& elements);
	void getElements(int ix, int iy, int iz, std::vector<T*>& elements);

	// Get grid infos
	unsigned int getNbCells() const { return _grid.size(); }
	unsigned int getResX() const { return _resX; }
	unsigned int getResY() const { return _resY; }
	unsigned int getResZ() const { return _resZ; }
	double getCellSize() const { return _cellSize; }
	Vector3DF getVolumeStart() const { return _volMin; }

	bool isCellEmpty(int ix, int iy, int iz) const { return _grid[getGridIndex(ix, iy, iz)].empty(); }

private:
	int getGridIndex(int ix, int iy, int iz) const;
	int getXIndex(const Vector3DF& position) const;
	int getYIndex(const Vector3DF& position) const;
	int getZIndex(const Vector3DF& position) const;
	int getXIndex(double xPos) const;
	int getYIndex(double yPos) const;
	int getZIndex(double zPos) const;

private:
	std::vector<std::vector<T> >	_grid;
	
	unsigned int	_resX;
	unsigned int	_resY;
	unsigned int	_resZ;
	double			_cellSize;
	Vector3DF		_volMin;
};	


// Include definitions
#include "SpatialGrid.hpp"

#endif	// SPATIALGRID_H
