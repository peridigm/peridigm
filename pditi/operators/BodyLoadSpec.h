/*
 * BodyLoadSpec.h
 *
 *  Created on: Jun 18, 2010
 *      Author: jamitch
 */

#ifndef BODYLOADSPEC_H_
#define BODYLOADSPEC_H_
#include "Loader.h"
#include "StageFunction.h"
#include <tr1/memory>

using std::tr1::shared_ptr;

namespace PdImp {
class BodyLoadSpec {
public:
	BodyLoadSpec(double unitDirection[3], const Pd_shared_ptr_Array<int>& pointIds);
	const Pd_shared_ptr_Array<int>& getPointIds() const { return localIds; }
	shared_ptr<PdImp::Loader> getStageLoader(const StageFunction& stageFunction) const;

private:
	/*
	 * Set of point Ids -- these Ids should be relevant to the fields passed in to the loader
	 */
	const Pd_shared_ptr_Array<int> localIds;

	/*
	 * Unit Direction
	 */
	double u[3];

	class BodyLoader;
	friend class BodyLoader;
	class BodyLoader : public PdImp::Loader {
	public:
		BodyLoader(const BodyLoadSpec& spec, const StageFunction& stageFunction);
		virtual void computeOwnedExternalForce(double lambda, Field_NS::Field<double> force) const;

	private:
		StageFunction stageFunction;
		const Pd_shared_ptr_Array<int> localIds;
		double u[3];

	};
};
}
#endif /* BODYLOADSPEC_H_ */
