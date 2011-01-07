/*
 * PimpStage.h
 *
 *  Created on: Jun 18, 2010
 *      Author: jamitch
 */

#ifndef PIMPSTAGE_H_
#define PIMPSTAGE_H_

#include <vector>
#include <tr1/memory>
using std::vector;
using std::tr1::shared_ptr;

namespace PdImp {

class Stage {
private:
	vector< shared_ptr<Loader> > loaders;

public:
	Stage() {}
	const vector< shared_ptr<Loader> >& getLoaders() const { return loaders; }
	void addLoader(shared_ptr<Loader>& loaderPtr);

};

}


#endif /* PIMPSTAGE_H_ */
