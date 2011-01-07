#ifndef KINEMATICBCSTAGE_H
#define KINEMATICBCSTAGE_H

namespace PdImp {

class StageFunction {

private:
	double start, end;

public:
	StageFunction() : start(0), end(0) {}

	StageFunction(double start, double end) : start(start), end(end){}

	/**
	 * Control function which produces step or proportional loading;
	 * Step loading is given with start==end
	 * Proportional loading is given by start != end
	 * @param lambda Load parameter for stage; it is assumed that 0<=lambda <=1.0
	 * @return Function value
	 */
	double value(double lambda) const {
		return (1-lambda)*start + lambda*end;
	}

	/**
	 * @return slope of function
	 */
	double slope() const { return end-start; }


	/**
	 * StageFunction for next stage
	 * @return new StageFunction will hold constant
	 */
	StageFunction next() const {
		double fStartNew = this->end;
		double fEndNew = fStartNew;
		return StageFunction(fStartNew,fEndNew);
	}



	/**
	 * StageFunction for next stage
	 * @param <code>endVal </code>value of loader at end of load step
	 * @return new proportional<code>StageFunction</code> with starting value of this end value and end value <code>endVal</code>
	 */
	StageFunction next(double endVal) const {
		double fStartNew = this->end;
		double fEndNew = endVal;
		return StageFunction(fStartNew,fEndNew);
	}

};

}
#endif
