#ifndef MESH_H
#define MESH_H

#include <Math/vector.h>
using Math::Vector3DF;

#include <vector>

class Mesh
{
public:
	struct Triangle
	{
		Triangle() {}
		Triangle(unsigned int v0, unsigned int v1, unsigned int v2)
		{
			v[0] = v0;
			v[1] = v1;
			v[2] = v2;
		}

		unsigned int v[3];
	};

public:
	Mesh() {}
	~Mesh() {}

	inline std::vector<Triangle>& triangles() { return _triangles; }
	inline std::vector<Vector3DF>& points() { return _points; }

	inline const std::vector<Triangle>& triangles() const  { return _triangles; }
	inline const std::vector<Vector3DF>& points() const { return _points; }


	inline void clear() { _triangles.clear(); _points.clear(); }
	inline void clearMemory()
	{ 
		std::vector<Triangle>().swap(_triangles);
		std::vector<Vector3DF>().swap(_points);
	}

private:
	std::vector<Triangle>	_triangles;
	std::vector<Vector3DF>	_points;
};

#endif	// MESH_H
