/**
 * Code de François Dagenais
 * Ce code est base sur le code de HoudiniSPH12 écrit par Jonathan Gagnon
 *
 */
#include <UT/UT_DSOVersion.h>

#include <GU/GU_PrimPoly.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_Director.h>

#include <GA/GA_AttributeRef.h>
#include <GEO/GEO_PrimPart.h>

#include "HoudiniSurfaceFromPoints.h"

// -----------------------------------------------------------------------------
// Create plugin inside houdini
// -----------------------------------------------------------------------------
void newSopOperator(OP_OperatorTable *table) {
     table->addOperator(new OP_Operator("pointsSurface",
					"MokkoSurfaceFromPoints",
					 HoudiniSurfaceFromPoints::myConstructor,
					 HoudiniSurfaceFromPoints::myTemplateList,
					 1,
					 1,
					 0));
}

// -----------------------------------------------------------------------------
// Plugin construction/destruction
// -----------------------------------------------------------------------------
static PRM_Name names[] = {
	PRM_Name("resolution", "MC cubes dimensions"),										// 0
	PRM_Name("influenceRadius", "Influence radius used by scalar function"),			// 1
	PRM_Name("particleRadius", "Desired distance of the surface from the particles"),	// 2
	PRM_Name("surfaceOnly", "Used near surface optimization"),							// 3
	PRM_Name("debugPoints", "Debug: Output points"),									// 4
	PRM_Name("debugSurfaceVertices", "Debug: Output surface vertices"),					// 5
	PRM_Name("useColorField", "Use color field calculations for near surface optimization (otherwise use cells)"),	// 6
	PRM_Name("computeNormals", "Compute normals"),										// 7
};

static PRM_Range positiveRange = PRMpositiveRange;

static PRM_Default defaultValues[] = {
	PRM_Default(0.1),										// 0
	PRM_Default(0.4),										// 1
	PRM_Default(0.1),										// 2
	PRM_Default(0),											// 3
	PRM_Default(0),											// 4
	PRM_Default(0),											// 5
	PRM_Default(0),											// 6
	PRM_Default(1)											// 7
};

PRM_Template HoudiniSurfaceFromPoints::myTemplateList[] = { 
	PRM_Template(PRM_FLT_J,	1, &names[0], &defaultValues[0], 0, &positiveRange),
	PRM_Template(PRM_FLT_J,	1, &names[1], &defaultValues[1], 0, &positiveRange),
	PRM_Template(PRM_FLT_J, 1, &names[2], &defaultValues[2], 0, &positiveRange),
	PRM_Template(PRM_TOGGLE, 1, &names[3], &defaultValues[3]),
	PRM_Template(PRM_TOGGLE, 1, &names[6], &defaultValues[6]),
	PRM_Template(PRM_TOGGLE, 1, &names[4], &defaultValues[4]),
	PRM_Template(PRM_TOGGLE, 1, &names[5], &defaultValues[5]),
	PRM_Template(PRM_TOGGLE, 1, &names[7], &defaultValues[7]),
	PRM_Template() };

OP_Node *
HoudiniSurfaceFromPoints::myConstructor(OP_Network *net, const char *name, OP_Operator *op) {
	return new HoudiniSurfaceFromPoints(net, name, op);
}

HoudiniSurfaceFromPoints::HoudiniSurfaceFromPoints(OP_Network *net, const char *name, OP_Operator *op) :
    SOP_Node(net, name, op),
	_hasSurfaceOnlyInformation(false)
{
}

HoudiniSurfaceFromPoints::~HoudiniSurfaceFromPoints() {
}

// -----------------------------------------------------------------------------
// Cooking! ( & stuff)
// -----------------------------------------------------------------------------
OP_ERROR HoudiniSurfaceFromPoints::cookMySop(OP_Context &context) {

    // Before we do anything, we must lock our inputs.  Before returning,
    //	we have to make sure that the inputs get unlocked.
	if (lockInputs(context) >= UT_ERROR_ABORT)
		return error();

	// Get Parameters
	fpreal now = context.getTime();
	double resolution = evalFloat("resolution", 0, now);
	double influenceRadius = evalFloat("influenceRadius", 0, now);
	double particleRadius = evalFloat("particleRadius", 0, now);
	bool surfaceOnly = (evalInt("surfaceOnly", 0, now) == 1);
	bool useColorField = (evalInt("useColorField", 0, now) == 1);
	bool computeNormals = (evalInt("computeNormals", 0, now) == 1);

	bool debugPoints = (evalInt("debugPoints", 0, now) == 1);
	bool debugSurfaceVertices = (evalInt("debugSurfaceVertices", 0, now) == 1);

	// Get input points
	std::vector<Vector3DF>	inputPoints;
	std::vector<double>		pointMasses;
	std::vector<double>		pointDensities;
	std::vector<double>		particlesSPHRadius;
	getPoints(context, inputPoints, pointMasses, pointDensities, particlesSPHRadius);
	if (surfaceOnly && useColorField && !_hasSurfaceOnlyInformation)
	{
		useColorField = false;

		this->addWarning(SOP_MESSAGE, "One or more of the necessary attributes (mass, density and/or pscale) are missing! Cannot use Color field to find surface points");
	}

	// Make sure have at least one point
	if (inputPoints.size()<1)
	{
		return error();
	}
	
	// Generate surface
	Mesh surface;
	std::vector<Vector3DF> normals;
	SurfaceFromPoints surfacer;
	if (surfaceOnly)
	{
		// Here, we assume every particle has the same sph radius
		double sphRadius = particlesSPHRadius[0];
		surfacer.createSurfaceFromPointsFaster(inputPoints,
										 	   surface,
											   normals,
										 	   resolution,
										 	   influenceRadius,
										 	   particleRadius,
											   pointMasses,
											   pointDensities,
											   sphRadius,
											   useColorField,
											   computeNormals);
	}
	else
	{
		surfacer.createSurfaceFromPoints(inputPoints,
										 surface,
										 normals,
										 resolution,
										 influenceRadius,
										 particleRadius,
										 computeNormals);
	}

	// Output geometry
	gdp->clearAndDestroy();	// Start from scratch
	if (debugPoints)
	{
		if (useColorField)
		{
			debugOutputPoints(inputPoints, surfacer._colorGradientNorm);
		}
		else
		{
			debugOutputPoints(inputPoints, surfacer._surfacePoints);
		}
	}
	else if (debugSurfaceVertices)
	{
		std::vector<double> dummy;
		debugOutputPoints(surfacer._surfaceVertices, dummy);
	}
	else
	{
		outputGeometry(surface, computeNormals, normals);
	}

	unlockInputs();
	return error();
}

const char *
HoudiniSurfaceFromPoints::inputLabel(unsigned) const
{
    return "Input geometry";
}

unsigned int HoudiniSurfaceFromPoints::disableParms()
{
	unsigned int changed = 0;
	float now = CHgetEvalTime();

	bool surfaceOnly = (evalInt("surfaceOnly", 0, now) == 1);
	changed = enableParm("useColorField", surfaceOnly);
	changed = enableParm("debugPoints", surfaceOnly);
	changed = enableParm("debugSurfacesVertices", surfaceOnly);

	return changed;
}

void HoudiniSurfaceFromPoints::getPoints(OP_Context&				context,
										 std::vector<Vector3DF>&	points,
										 std::vector<double>&		masses,
										 std::vector<double>&		densities,
										 std::vector<double>&		sphRadius)
{
	// Get source
	const GU_Detail* source = inputGeo(0, context);
	GEO_PointList pnts = source->points();

	GA_ROAttributeRef srcMassAttr = source->findFloatTuple(GA_ATTRIB_POINT, "mass", 1);
	GA_ROAttributeRef srcDensityAttr = source->findFloatTuple(GA_ATTRIB_POINT, "density", 1);
	GA_ROAttributeRef srcSPHRadius = source->findFloatTuple(GA_ATTRIB_POINT, "pscale", 1);

	_hasSurfaceOnlyInformation = (srcMassAttr.isValid() && srcDensityAttr.isValid()
									&& srcSPHRadius.isValid());

	// Allocate memory for points and properties
	points.resize(pnts.entries());
	if (_hasSurfaceOnlyInformation)
	{
		masses.resize(pnts.entries());
		densities.resize(pnts.entries());
		sphRadius.resize(pnts.entries());
	}

	// Get points
	GEO_Point ppt(source->getPointMap(), GA_INVALID_OFFSET);
	Vector3DF point;
	points.resize(pnts.entries());
	for (int p=0; p<pnts.entries(); ++p)
	{
		// TODO: Should use source->getGEOPoint!!!
		ppt = GEO_Point(source->getPointMap(), source->pointOffset(p));

		point.x = ppt.getPos3().x();
		point.y = ppt.getPos3().y();
		point.z = ppt.getPos3().z();

		points[p] = point;

		if (_hasSurfaceOnlyInformation)
		{
			masses[p] = ppt.getValue<float>(srcMassAttr);
			densities[p] = ppt.getValue<float>(srcDensityAttr);	
			sphRadius[p] = ppt.getValue<float>(srcSPHRadius);
		}
	}
}

void HoudiniSurfaceFromPoints::outputGeometry(const Mesh&					mesh,
											  bool							outputNormals,
											  const std::vector<Vector3DF>&	normals)
{
	const std::vector<Vector3DF>& points = mesh.points();
	const std::vector<Mesh::Triangle>& triangles = mesh.triangles();

	GA_RWAttributeRef normalAttr;
	if (outputNormals)
	{
		normalAttr = gdp->addFloatTuple(GA_ATTRIB_POINT, "N", 3);
	}

	// Generate points
	for (int p=0; p<points.size(); ++p)
	{
		const Vector3DF& point = points[p];

		GEO_Point *houdiniPoint = gdp->appendPointElement();
		houdiniPoint->setPos3(UT_Vector3(point.x, point.y, point.z));

		if (outputNormals)
		{
			const Vector3DF& n = normals[p];

			houdiniPoint->setValue<UT_Vector3>(normalAttr, UT_Vector3(n.x, n.y, n.z));
		}
	}

	// Generate triangles
	for (int t=0; t<triangles.size(); ++t)
	{
		// Create primitive
		GEO_PrimPoly *prim = dynamic_cast<GEO_PrimPoly*>(gdp->appendPrimitive(GEO_PRIMPOLY));
		prim->setSize(0);

		// Generate vertices
		prim->appendVertex(gdp->points()[triangles[t].v[0]]);
		prim->appendVertex(gdp->points()[triangles[t].v[2]]);
		prim->appendVertex(gdp->points()[triangles[t].v[1]]);

		// End primitive (this will generate a polygone instead of a polyline)
		prim->close();
	}
}

void HoudiniSurfaceFromPoints::debugOutputPoints(const std::vector<Vector3DF>&	points,
												 const std::vector<double>&		values)
{
	GEO_PrimParticle *mySystem;
	GEO_Point ppt(gdp->getPointMap(), GA_INVALID_OFFSET);
	GA_RWAttributeRef	outputAttr;

	bool outputValues = (values.size()>=points.size());
	if (outputValues)
	{
		outputAttr = gdp->addFloatTuple(GA_ATTRIB_POINT, "colorGradient", 1);
	}

	mySystem = (GEO_PrimParticle *) gdp->appendPrimitive(GEO_PRIMPART);
	mySystem->clearAndDestroy();

	for (int p=0; p<points.size(); ++p)
	{
		GA_Offset vtxoff = mySystem->giveBirth();

		ppt= GEO_Point(gdp->getPointMap(), gdp->vertexPoint(vtxoff));

		UT_Vector3 pos;
		pos.x() = points[p].x;
		pos.y() = points[p].y;
		pos.z() = points[p].z;

		ppt.setPos(pos);

		if (outputValues)
		{
			ppt.setValue<float>(outputAttr, values[p]);
		}
	}
}
