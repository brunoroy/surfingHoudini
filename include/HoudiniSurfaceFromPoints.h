#ifndef HOUDINISURFACEFROMPOINTS_H
#define HOUDINISURFACEFROMPOINTS_H

#include "SurfaceFromPoints.h"

#include <Math/vector.h>
using Math::Vector3DF;

#include <SOP/SOP_Node.h>
#include <memory>
#include <vector>

using namespace std;

class HoudiniSurfaceFromPoints : public SOP_Node {
public:
	HoudiniSurfaceFromPoints(OP_Network *net, const char *name, OP_Operator *op);
	virtual ~HoudiniSurfaceFromPoints();

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

	void outputGeometry(const							Mesh& mesh,
						bool							outputNormals,
						const std::vector<Vector3DF>&	normals);

	void debugOutputPoints(const std::vector<Vector3DF>&	points,
						   const std::vector<double>&		values);

	bool _hasSurfaceOnlyInformation;
};

#endif	// HOUDINISURFACEFROMPOINTS_H
