#ifndef __cache_object__
#define __cache_object__

#include <functional>

template <typename T>
class cache_object {
public:
	cache_object(const std::function < T (void)>& f)
		: _isGood(false), _f(f)
	{
	}

	cache_object(void)
		: _isGood(false)
	{}

	T value() const {
		if (!_isGood) {
			_object = _f();
			_isGood = true;
		}
		return _object;
	}
private:
	mutable T _object;
	mutable bool _isGood;
	const std::function <T(void)> _f;
};

template<typename T>
cache_object<T> make_cache_func(const std::function <T(void)>& f)
{
	return cache_object<T>(f);
}

#endif
