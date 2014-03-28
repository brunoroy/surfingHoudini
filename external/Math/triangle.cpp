
#include "triangle.h"

using namespace Math;

bool SameSide(Vector3DF p1,Vector3DF p2,Vector3DF a, Vector3DF b)
{
	/*
	UT_Vector3 cp1 = cross(b-a, p1-a);
	UT_Vector3 cp2 = cross(b-a, p2-a);
	if (dot(cp1, cp2) >= 0)
	 return true;
	else return false;
	*/

	Vector3DF ab = b - a;
	Vector3DF ap1 = p1 -a;
	Vector3DF ap2 = p2 - a;
	Vector3DF cp1 = ab.Cross(ap1);
    	Vector3DF cp2 = ab.Cross(ap2);
    	if (cp1.Dot(cp2) >= 0)
		return true;
   	else return false;
}

float TriangleAir(Vector3DF a,Vector3DF b,Vector3DF c)
{
	/*
	UT_Vector3 x =(c-b);
	UT_Vector3 y =(c-a);
	UT_Vector3 z = cross(x,y);
	
	return 0.5* z.length();
	*/

	Vector3DF x =(c-b);
	Vector3DF y =(c-a);
	Vector3DF z = x.Cross(y);
	
	return 0.5* z.Length();
}


Vector3DF sous(Vector3DF A, Vector3DF B)
{
	Vector3DF r;
	r.x = B.x - A.x;
	r.y = B.y - A.y;
	r.z = B.z - A.z;
	return r;
}

Vector3DF add(Vector3DF A, Vector3DF B)
{
	Vector3DF r;
	r.x = B.x + A.x;
	r.y = B.y + A.y;
	r.z = B.z + A.z;
	return r;
}


bool Triangle::PointInTriangle2d(Vector3DF p, Vector3DF a,Vector3DF b,Vector3DF c)
{
    if (SameSide(p,a, b,c) && SameSide(p,b, a,c) && SameSide(p,c, a,b))
	return true;
    else return false;
}

Vector3DF Triangle::CreateVector( Vector3DF a,Vector3DF b)
{
	Vector3DF result;
	result.x = b.x - a.x;
	result.y = b.y - a.y;
	result.z = b.z - a.z;
	return result;
}

float Triangle::TriangleAir(Vector3DF a,Vector3DF b,Vector3DF c)
{
	Vector3DF x =(c-b);
	Vector3DF y =(c-a);
	Vector3DF z = x.Cross(y);
	
	return 0.5* z.Length();
}

bool Triangle::IsPointInTriangle2D(const Vector3DF& p, const Vector3DF& a, const Vector3DF& b, const Vector3DF& c)
{
	// Determine on wich side of each lines the point lies
	double sideAB = (p.x-a.x)*(b.y-a.y) - (p.y-a.y)*(b.x-a.x);
	double sideBC = (p.x-b.x)*(c.y-b.y) - (p.y-b.y)*(c.x-b.x);
	double sideCA = (p.x-c.x)*(a.y-c.y) - (p.y-c.y)*(a.x-c.x);

	// Are we always on the good side of the line?
	return (sideAB>=0.0 && sideBC>=0.0 && sideCA>=0.0);
}

bool Triangle::IsPointInTriangle(Vector3DF  p, Vector3DF a,Vector3DF b,Vector3DF c)
{
	/*
		// Compute vectors
		UT_Vector3 v0 = t->c - t->a;
		UT_Vector3 v1 = t->b - t->a;
		UT_Vector3 v2 = p - t->a;

		// Compute dot products
		float dot00 = dot(v0, v0);
		float dot01 = dot(v0, v1);
		float dot02 = dot(v0, v2);
		float dot11 = dot(v1, v1);
		float dot12 = dot(v1, v2);

		// Compute barycentric coordinates
		float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
		float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
		float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

		// Check if point is in triangle
		return (u >= 0) && (v >= 0) && (u + v <= 1);
	*/
	
	float epsilon = 0.00001;
	Vector3DF v0 = c - a;
	Vector3DF v1 = b - a;
	Vector3DF v2 = p - a;


	// Compute dot products
	float dot00 = v0.Dot(v0);
	float dot01 = v0.Dot(v1);
	float dot02 = v0.Dot(v2);
	float dot11 = v1.Dot(v1);
	float dot12 = v1.Dot(v2);

	// Compute barycentric coordinates
	float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point is in triangle
	return (u >= 0) && (v >= 0) && (u + v <= 1);
}

bool Triangle::IsPointInTriangle2(Vector3DF  p, Vector3DF a,Vector3DF b,Vector3DF c)
{
	/*
		// Compute vectors
		UT_Vector3 v0 = t->c - t->a;
		UT_Vector3 v1 = t->b - t->a;
		UT_Vector3 v2 = p - t->a;

		// Compute dot products
		float dot00 = dot(v0, v0);
		float dot01 = dot(v0, v1);
		float dot02 = dot(v0, v2);
		float dot11 = dot(v1, v1);
		float dot12 = dot(v1, v2);

		// Compute barycentric coordinates
		float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
		float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
		float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

		// Check if point is in triangle
		return (u >= 0) && (v >= 0) && (u + v <= 1);
	*/

	float epsilon = 0.0001;
	Vector3DF v0 = c - a;
	Vector3DF v1 = b - a;
	Vector3DF v2 = p - a;


	Vector3DF e1 = b-a;
	Vector3DF e2 = b-c;
	Vector3DF e3 = a-c;

	Vector3DF n1 = e1;
	Vector3DF n2 = e2;
	Vector3DF n3 = e3;

	n1.Normalize();
	n2.Normalize();
	n3.Normalize();

	double d1 = squaredDistToEdge(p, a, n1, e1.Length() );
	double d2 = squaredDistToEdge(p, b, n2, e2.Length() );
	double d3 = squaredDistToEdge(p, c, n3, e3.Length() );


	// Compute dot products
	float dot00 = v0.Dot(v0);
	float dot01 = v0.Dot(v1);
	float dot02 = v0.Dot(v2);
	float dot11 = v1.Dot(v1);
	float dot12 = v1.Dot(v2);

	// Compute barycentric coordinates
	float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	// Check if point is in triangle
	return (((u >= 0) && (v >= 0) && (u + v <= 1)) || d1 < epsilon || d2 < epsilon || d3 < epsilon);
	//return (u >= 0) && (v >= 0) && (u + v <= 1);
}
bool Triangle::IsPointInTriangleRange(Vector3DF  p, Vector3DF a,Vector3DF b,Vector3DF c, float s)
{

	//-------------- scaling ------------Triangle::----
	
	
	Vector3DF center = (a+b+c);//3;
	center /= 3;
	Vector3DF s1 = a - center;
	Vector3DF s2 = b - center;
	Vector3DF s3 = c - center;

	//s1.normalize();
	//s2.normalize();
	//s3.normalize();

	s1 *=s;
	s2 *=s;
	s3 *=s;		

	a = a+(s1);
	b = b+(s2);
	c = c+(s3);
	
	//----------------------------------------

	
  	Vector3DF l1,l2,l3;
	float angle;
	float dot1,dot2,dot3;
	l1 = p- a;
	l2 = p- b;
	l3 = p- c;
	l1.Normalize();
	l2.Normalize();
	l3.Normalize();
	
	dot1 = l1.Dot(l2);
	dot2 = l2.Dot(l3);
	dot3 = l3.Dot(l1);

	angle = acos(dot1) + acos(dot2) +acos(dot3);

	if ( (angle < ((2*3.14159)+0.01)) && (angle > ((2*3.14159)-0.01)) )
		return true;

	return false;

	
}

Vector3DF Triangle::GetBarycentricPosition(Vector3DF A,Vector3DF B, Vector3DF C, Vector3DF a, Vector3DF b, Vector3DF c, Vector3DF position)
{

	float totalAir = TriangleAir(A,B,C);
	float apbAir = TriangleAir(A,B,position);
	float bpcAir = TriangleAir(B,C,position);
	float apcAir = TriangleAir(A,C, position);

	float alpha = bpcAir/totalAir;
	float beta =  apcAir/totalAir;
	float gamma = apbAir/totalAir;

	a *= alpha;
	b *= beta;
	c *= gamma;

	return a+b+c;


}

Vector3DF Triangle::ProjectPointOnTriangle(Vector3DF p,  Vector3DF a,Vector3DF b,Vector3DF c)
{
	/*

	UT_Vector3 v1 = a - b;
	UT_Vector3 v2 = a - c;
	UT_Vector3 n = cross(v1,v2);GRID
	n.normalize();	

	UT_Vector3 r0 = a;
	UT_Vector3 v3 = r0 - p;
	UT_Vector3 v3m = v3-n*dot(n,v3);
	UT_Vector3 rm = r0-v3m;
	return rm;
	*/
	Vector3DF v1 = b - a;
	Vector3DF v2 = c - a;
	Vector3DF n = v1.Cross(v2);
	n.Normalize();	

	Vector3DF r0 = a;
	Vector3DF v3 = p - r0;
	float dot = n.Dot(v3);
	Vector3DF ndot = n;
	ndot *= dot;
	Vector3DF v3m = v3-ndot;
	Vector3DF rm = r0 + v3m;
	return rm;
}


Vector3DF Triangle::GetUVBarycenter(Vector3DF a,Vector3DF b, Vector3DF c,Vector3DF A,Vector3DF B, Vector3DF C, Vector3DF P)
{
	/*    
	UT_Vector3 v0 = c - a;
	UT_Vector3 v1 = b - a;
	UT_Vector3 v2 = P - a;

	// Compute dot products
	float dot00 = dot(v0, v0);
	float dot01 = dot(v0, v1);
	float dot02 = dot(v0, v2);
	float dot11 = dot(v1, v1);
	float dot12 = dot(v1, v2);

	// Compute barycentric coordinatesTriangle::
	float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	return A + u*(C-A) + v* (B-A);
	*/
	
	// Compute vectors        
	//Vector3DF v0 = c - a;
	//Vector3DF v1 = b - a;
	//Vector3DF v2 = P - a;

	Vector3DF v0 = sous(a,c); 
	Vector3DF v1 = sous(a,b);
	Vector3DF v2 = sous(a,P);


	
	// Compute dot products
	float dot00 = v0.Dot(v0);
	float dot01 = v0.Dot(v1);
	float dot02 = v0.Dot(v2);
	float dot11 = v1.Dot(v1);
	float dot12 = v1.Dot(v2);

	// Compute barycentric coordinates
	float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
	float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	//return A + u*(C-A) + v* (B-A);
	Vector3DF ac = sous(A,C);
	Vector3DF ab = sous(A,B);
	ac *= u;
	ab *= v;
	Vector3DF r;
	return A + ac + ab;
}

void Triangle::BoundingBox(Vector3DF a,Vector3DF b,Vector3DF c, Vector3DF &min,Vector3DF &max)
{
	min.x = a.x;
	min.y = a.y;
	min.z = a.z;
	max = min;

	// check min

	if (min.x > b.x)
	{
		min.x = b.x;
	}
	if (min.y > b.y)
	{
		min.y = b.y;
	}
	if (min.z > b.z)
	{
		min.z = b.z;
	}

	if (min.x > c.x)
	{
		min.x = c.x;
	}
	if (min.y > c.y)
	{
		min.y = c.y;
	}
	if (min.z > c.z)
	{
		min.z = c.z;
	}

	//check max

	if (max.x < b.x)
	{
		max.x = b.x;
	}
	if (max.y < b.y)
	{
		max.y = b.y;
	}
	if (max.z < b.z)
	{
		max.z = b.z;
	}

	if (max.x < c.x)
	{
		max.x = c.x;
	}
	if (max.y < c.y)
	{
		max.y = c.y;
	}
	if (max.z < c.z)
	{
		max.z = c.z;
	}
}

