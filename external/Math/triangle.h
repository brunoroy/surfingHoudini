#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include "vector.h"


using namespace std;


namespace Math {

class Triangle{
	public :

	bool static SameSide(Vector3DF p1,Vector3DF p2,Vector3DF a, Vector3DF b);
	bool static PointInTriangle2d(Vector3DF p, Vector3DF a,Vector3DF b,Vector3DF c);
	Vector3DF static CreateVector( Vector3DF a,Vector3DF b);
	float static TriangleAir(Vector3DF a,Vector3DF b,Vector3DF c);
	bool static IsPointInTriangle(Vector3DF  p, Vector3DF a,Vector3DF b,Vector3DF c);
	bool static IsPointInTriangle2(Vector3DF  p, Vector3DF a,Vector3DF b,Vector3DF c);
	bool static IsPointInTriangle2D(const Vector3DF& p, const Vector3DF& a, const Vector3DF& b, const Vector3DF& c);
	bool static IsPointInTriangleRange(Vector3DF  p, Vector3DF a,Vector3DF b,Vector3DF c, float s);
	Vector3DF static ProjectPointOnTriangle(Vector3DF p,  Vector3DF a,Vector3DF b,Vector3DF c);
	Vector3DF static GetUVBarycenter(Vector3DF a,Vector3DF b, Vector3DF c,Vector3DF A,Vector3DF B, Vector3DF C, Vector3DF P);
	Vector3DF GetBarycentricPosition(Vector3DF A,Vector3DF B, Vector3DF C, Vector3DF a, Vector3DF b, Vector3DF c, Vector3DF position);
	inline
		static double squaredDistToEdge(const Vector3DF& p,
								 const Vector3DF& v0,
								 const Vector3DF& edgeDirection,
								 double			  edgeLength)
		{
			// Project point onto the edge
			Vector3DF v0p(p);
			v0p -= v0;

			double d = v0p.x*edgeDirection.x + v0p.y*edgeDirection.y + v0p.z*edgeDirection.z;

			if (d<0.0) d = 0.0;
			if (d>edgeLength) d = edgeLength;

			Vector3DF projP(edgeDirection);
			projP *= d;

			// Compute the squared distance
			Vector3DF delta(projP);
			delta -= v0p;

			return delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
		}
	void static BoundingBox(Vector3DF a,Vector3DF b,Vector3DF c, Vector3DF &min,Vector3DF &max);
};

}

#endif
