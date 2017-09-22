#ifndef __doubleError_H__
#define __doubleError_H__

#include <cmath>
#include <ostream>

// Object to try to encapsulate error propagation
class doubleError {
public:
	doubleError(double centralValue = 0, double error = 0.0)
		: _v(centralValue), _err2(error*error)
	{}
	// Internal, to avoid squaring.
	doubleError(bool doit, double value, double err2)
		: _v(value), _err2(err2)
	{}


	double value() const { return _v; }
	double err() const { return sqrt(_err2); }
	double err2() const { return _err2; }
	double ferr() const {
		if (_v == 0)
			return 0;
		return err() / _v;
	}

	// Add constant, update in place.
	doubleError &operator+= (const doubleError &e2)
	{
		_v += e2._v;
		_err2 += e2._err2;
		return *this;
	}

	// Add offset - update in place.
	doubleError &operator+= (const double e2)
	{
		_v += e2;
		return *this;
	}

	// Multiply
	doubleError &operator*= (const doubleError &e2)
	{
		double newv = _v * e2.value();
		double ferr1 = ferr();
		double ferr2 = e2.ferr();
		double ferrSquared = ferr1*ferr1 + ferr2*ferr2;
		_v = newv;
		_err2 = newv*newv * ferrSquared;

		return *this;
	}

	// Simple equality test.
	bool operator== (const doubleError &other) {
		return _v == other._v
			&& _err2 == other._err2;
	}

	// And the opposite
	bool operator!= (const doubleError &other) {
		return !(*this == other);
	}

private:
	double _v;
	double _err2;

	doubleError scaleBy(double factor)
	{
		return doubleError(_v * factor, _err2 * factor);
	}
};

std::ostream &operator<< (std::ostream &s, const doubleError &e)
{
	s << e.value() << " +- " << e.err();
	return s;
}

doubleError operator/ (double e1, const doubleError &e2)
{
	double newv = e1 / e2.value();
	return doubleError(newv, newv * e2.ferr());
}

// Add constant
doubleError operator+ (const doubleError &e1, const doubleError &e2)
{
	return doubleError(false, e1.value() + e2.value(), e1.err2() + e2.err2());
}

// Add offset
doubleError operator+(const doubleError &e1, double e2)
{
	return doubleError(false, e1.value() + e2, e1.err2());
}

// Multiply
doubleError operator*(const doubleError &e1, double e2)
{
	double newv = e1.value() * e2;
	return doubleError(newv, newv * e1.ferr());
}

doubleError operator*(const doubleError &e1, const doubleError &e2)
{
	double newv = e1.value() * e2.value();
	double ferr1 = e1.ferr();
	double ferr2 = e2.ferr();
	double ferrSquared = ferr1*ferr1 + ferr2*ferr2;
	return doubleError(newv, newv * sqrt(ferrSquared));
}

// Divide offset
doubleError operator/ (const doubleError &e1, double e2)
{
	return e1 * (1.0 / e2);
}

// Real divide
doubleError operator/ (const doubleError &e1, const doubleError &e2)
{
	return e1 * (1.0 / e2);
}

#endif
