//#include <UT/UT_DSOVersion.h>

#include <GU/GU_PrimPoly.h>
#include <PRM/PRM_Include.h>
#include <OP/OP_Director.h>

#include <GA/GA_AttributeRef.h>
#include <GEO/GEO_PrimPart.h>

#include "surfingHoudini.h"

// -----------------------------------------------------------------------------
// Create plugin inside houdini
// -----------------------------------------------------------------------------
void newSopOperator(OP_OperatorTable *table) {
     table->addOperator(new OP_Operator("surfingSAT",
                    "SurfingSAT",
                     SurfingHoudini::myConstructor,
                     SurfingHoudini::myTemplateList,
					 1,
					 1,
					 0));
}

// -----------------------------------------------------------------------------
// Plugin construction/destruction
// -----------------------------------------------------------------------------
static PRM_Name names[] = {
    PRM_Name("resolution", "Resolution")
};

static PRM_Range unitRange(PRM_RANGE_RESTRICTED, 0.01f, PRM_RANGE_RESTRICTED, 0.5f);
static PRM_Range positiveRange = PRMpositiveRange;

static PRM_Default defaultValues[] = {
    PRM_Default(0.1)
};

PRM_Template SurfingHoudini::myTemplateList[] = {
    PRM_Template(PRM_FLT_J,	1, &names[0], &defaultValues[0], 0, &unitRange),
    PRM_Template()
};

OP_Node *
SurfingHoudini::myConstructor(OP_Network *net, const char *name, OP_Operator *op) {
    return new SurfingHoudini(net, name, op);
}

SurfingHoudini::SurfingHoudini(OP_Network *net, const char *name, OP_Operator *op) :
    SOP_Node(net, name, op),
	_hasSurfaceOnlyInformation(false)
{
}

SurfingHoudini::~SurfingHoudini() {
}

// -----------------------------------------------------------------------------
// Cooking! ( & stuff)
// -----------------------------------------------------------------------------
OP_ERROR SurfingHoudini::cookMySop(OP_Context &context) {

    // Before we do anything, we must lock our inputs.  Before returning,
    //	we have to make sure that the inputs get unlocked.
	if (lockInputs(context) >= UT_ERROR_ABORT)
		return error();

	// Get Parameters
	fpreal now = context.getTime();
	double resolution = evalFloat("resolution", 0, now);

	// Get input points
	std::vector<Vector3DF>	inputPoints;
	std::vector<double>		pointMasses;
	std::vector<double>		pointDensities;
	std::vector<double>		particlesSPHRadius;
	getPoints(context, inputPoints, pointMasses, pointDensities, particlesSPHRadius);

	// Make sure have at least one point
	if (inputPoints.size()<1)
	{
		return error();
	}
	
	// Generate surface
	Mesh surface;
	std::vector<Vector3DF> normals;
	SurfaceFromPoints surfacer;

    surfacer.createSurfaceFromPoints(inputPoints, surface, resolution);

	// Output geometry
	gdp->clearAndDestroy();	// Start from scratch
    outputGeometry(surface);

	unlockInputs();
	return error();
}

const char *
SurfingHoudini::inputLabel(unsigned) const
{
    return "Input geometry";
}

unsigned int SurfingHoudini::disableParms()
{
	unsigned int changed = 0;

	return changed;
}

void SurfingHoudini::getPoints(OP_Context&				context,
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

void SurfingHoudini::outputGeometry(const Mesh& mesh)
{
	const std::vector<Vector3DF>& points = mesh.points();
	const std::vector<Mesh::Triangle>& triangles = mesh.triangles();

	// Generate points
	for (int p=0; p<points.size(); ++p)
	{
		const Vector3DF& point = points[p];

		GEO_Point *houdiniPoint = gdp->appendPointElement();
		houdiniPoint->setPos3(UT_Vector3(point.x, point.y, point.z));
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
