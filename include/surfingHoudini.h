#ifndef HOUDINI_SURFING_H
#define HOUDINI_SURFING_H

#include "SurfaceFromPoints.h"

#include <Math/vector.h>
using Math::Vector3DF;

#include <SOP/SOP_Node.h>
#include <memory>
#include <vector>

using namespace std;

class SurfingHoudini : public SOP_Node {
public:
    SurfingHoudini(OP_Network *net, const char *name, OP_Operator *op);
    virtual ~SurfingHoudini();

	static PRM_Template myTemplateList[];
	static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

	unsigned int disableParms();

protected:
    virtual const char* inputLabel(unsigned idx) const;
	virtual OP_ERROR cookMySop(OP_Context &context);

private:
	void getPoints(OP_Context&				context,
				   std::vector<Vector3DF>&	points,
				   std::vector<double>&		masses,
				   std::vector<double>&		densities,
				   std::vector<double>&		sphRadius);

    void outputGeometry(const Mesh& mesh);

	bool _hasSurfaceOnlyInformation;
};

#endif	// HOUDINI_SURFING_H
