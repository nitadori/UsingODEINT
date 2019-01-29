#include <boost/numeric/odeint.hpp>
#include <cstdio>
#include <cmath>

struct Vec2D{
	double x, y;

	double abs() const {
		return std::sqrt(*this * *this);
	}
	double operator*(const Vec2D &rhs) const {
		return x*rhs.x + y*rhs.y;
	}
	Vec2D operator+(const Vec2D &rhs) const {
		return {x+rhs.x, y+rhs.y};
	}
	Vec2D operator-(const Vec2D &rhs) const {
		return {x-rhs.x, y-rhs.y};
	}
	Vec2D &operator+=(const Vec2D &rhs) {
		return (*this = *this + rhs);
	}
	Vec2D &operator-=(const Vec2D &rhs) {
		return (*this = *this - rhs);
	}
	friend Vec2D operator*(const double s, const Vec2D &v){
		return {s*v.x, s*v.y};
	}
};

struct Threebody{
	static double mass[3];
	Vec2D  pos [3];
	Vec2D  vel [3];

	void operator()(const Threebody &val, Threebody &der, double t){
		auto grav = [](Vec2D ri, Vec2D rj)->Vec2D {
			Vec2D  rij = rj - ri;
			double r2 = rij * rij;
			double r2inv = 1.0 / r2;
			double r3inv = r2inv * std::sqrt(r2inv);

			return r3inv * rij;
		};
		Vec2D f01 = grav(val.pos[0], val.pos[1]);
		Vec2D f12 = grav(val.pos[1], val.pos[2]);
		Vec2D f20 = grav(val.pos[2], val.pos[0]);
		Vec2D f0 = mass[1] * f01 - mass[2] * f20;
		Vec2D f1 = mass[2] * f12 - mass[0] * f01;
		Vec2D f2 = mass[0] * f20 - mass[1] * f12;

		der.vel[0] = f0;
		der.vel[1] = f1;
		der.vel[2] = f2;

		der.pos[0] = val.vel[0];
		der.pos[1] = val.vel[1];
		der.pos[2] = val.vel[2];
	}

	double *begin(){
		return (double*)&pos;
	}
	double *end(){
		return begin() + 6;
	}
	const double *begin() const {
		return (double*)&pos;
	}
	const double *end() const {
		return begin() + 6;
	}
	void resize(const size_t) {}
};
double Threebody::mass[3];

namespace boost{ namespace numeric{ namespace odeint{
	template<>
	struct is_resizeable<Threebody>{
			static const bool value = boost::true_type::value;
	};
} } }

