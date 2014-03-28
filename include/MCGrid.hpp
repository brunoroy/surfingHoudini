struct MCGrid::MCVertex
{
	MCVertex() : value(1000000000000000000.0) { points[0] = -1; points[1] = -1; points[2] = -1; }

	double value;
	int points[3];
	int gridIndex;
	Vector3DF normal;
};

inline
void MCGrid::getCellsInRadius(const Vector3DF&	position,
							  double			radius,
							  unsigned int&		minX,
							  unsigned int&		maxX,
							  unsigned int&		minY,
							  unsigned int&		maxY,
							  unsigned int&		minZ,
							  unsigned int&		maxZ) const
{
	int xmin, ymin, zmin;	// Because unsigned int will loop if we go below zero

    xmin = static_cast<int>(floor( ((position.x-_volMin.x)-radius)/_cubeSize  ));
    ymin = static_cast<int>(floor( ((position.y-_volMin.y)-radius)/_cubeSize  ));
    zmin = static_cast<int>(floor( ((position.z-_volMin.z)-radius)/_cubeSize  ));

    maxX = static_cast<int>(floor( ((position.x-_volMin.x)+radius)/_cubeSize  ));
    maxY = static_cast<int>(floor( ((position.y-_volMin.y)+radius)/_cubeSize  ));
    maxZ = static_cast<int>(floor( ((position.z-_volMin.z)+radius)/_cubeSize  ));

	minX = (xmin<0) ? 0 : static_cast<unsigned int>(xmin);
	if (maxX >= _resX) maxX = _resX-1;

	minY = (ymin<0) ? 0 : static_cast<unsigned int>(ymin);
	if (maxY >= _resY) maxY = _resY-1;

	minZ = (zmin<0) ? 0 : static_cast<unsigned int>(zmin);
	if (maxZ >= _resZ) maxZ = _resZ-1;
}

inline
void MCGrid::getVertexPosition(unsigned int ix,
									  unsigned int iy,
									  unsigned int iz,
									  Vector3DF& position) const
{
	position = Vector3DF(ix, iy, iz);
	position *= _cubeSize;
    position += _volMin;
}

inline
double MCGrid::getScalarValue(unsigned int ix, unsigned int iy, unsigned int iz) const
{
	int dataIndex = _vertices[getGridIndex(ix, iy, iz)];
	return (dataIndex > -1) ? _verticesData[dataIndex].value : 1000000000000000000.0;
}

inline
void MCGrid::setScalarValue(unsigned int ix, unsigned int iy, unsigned int iz, double value)
{
	int gridIndex = getGridIndex(ix, iy, iz);
	int dataIndex = _vertices[gridIndex];

	// Create entry in data array if necessary
    if (dataIndex < 0)
	{
		dataIndex = _verticesData.size();
		_verticesData.push_back(MCVertex());
		_vertices[gridIndex] = dataIndex;

		_verticesData[dataIndex].gridIndex = gridIndex;
    }

	// Set value
	_verticesData[dataIndex].value = value;
}

inline
unsigned int MCGrid::getGridIndex(unsigned int ix, unsigned int iy, unsigned int iz) const
{
	return ix + (iy*_resX) + (iz*_resX*_resY);
}

inline unsigned int MCGrid::getXIndex(unsigned int gridIndex) const
{
	return gridIndex = gridIndex % _resX;
}

inline unsigned int MCGrid::getYIndex(unsigned int gridIndex) const
{
	return gridIndex = (gridIndex % (_resX*_resY)) / _resX;
}

inline unsigned int MCGrid::getZIndex(unsigned int gridIndex) const
{
	return gridIndex = gridIndex / (_resX*_resY);
}

inline
unsigned int MCGrid::getEdgePoint(MCVertex& v1,
								  MCVertex& v2,
								  MCVertex& v3,
								  MCVertex& v4,
								  MCVertex& v5,
								  MCVertex& v6,
								  MCVertex& v7,
								  MCVertex& v8,
								  int edgeNo,
								  std::vector<Vector3DF>& points) const
{
	MCVertex* va = 0x0;
	MCVertex* vb = 0x0;
	int axis = 0;

	// Get vertex and axis associated with the edge
    if (edgeNo == 1)
    {
		va = &v1;
		vb = &v2;
		axis = 0;
    }
	else if (edgeNo == 2)
	{
		va = &v2;
		vb = &v3;
		axis = 1;
	}
	else if (edgeNo == 3)
	{
		va = &v4;
		vb = &v3;
		axis = 0;
	}
    else if (edgeNo == 4)
	{
		va = &v1;
		vb = &v4;
		axis = 1;
	}
    else if (edgeNo == 5)
	{
		va = &v5;
		vb = &v6;
		axis = 0;
	}
	else if (edgeNo == 6)
	{
		va = &v6;
		vb = &v7;
		axis = 1;
	}
	else if (edgeNo == 7)
	{
		va = &v8;
		vb = &v7;
		axis = 0;
	}
	else if (edgeNo == 8)
	{
		va = &v5;
		vb = &v8;
		axis = 1;
	}
	else if (edgeNo == 9)
	{
		va = &v1;
		vb = &v5;
		axis = 2;
	}
	else if (edgeNo == 10)
	{
		va = &v2;
		vb = &v6;
		axis = 2;
	}
	else if (edgeNo == 11)
	{
		va = &v4;
		vb = &v8;
		axis = 2;
	}
	else //if (edgeNo == 12)
	{
		va = &v3;
		vb = &v7;
		axis = 2;
    }

	// Get edge point index
	int pointID = va->points[axis];
	if (pointID == -1)
	{
		// Create point if it hasn't been created yet
		// First compute position using linear interpolation
		// between the two vertices
		Vector3DF posA, posB, pos;
		getVertexPosition(getXIndex(va->gridIndex),
						  getYIndex(va->gridIndex),
						  getZIndex(va->gridIndex),
						  posA);
        getVertexPosition(getXIndex(vb->gridIndex),
                          getYIndex(vb->gridIndex),
                          getZIndex(vb->gridIndex),
						  posB);
		
		double t = (0.0 - va->value) / (vb->value - va->value);
        posA *= 1.0-t;
        posB *= t;

		pos = posA;
        pos += posB;

		// add point
		pointID = points.size();
		points.push_back(pos);
		va->points[axis] = pointID;

        //std::clog << points.size() << std::endl;
	}

	return pointID;
}

