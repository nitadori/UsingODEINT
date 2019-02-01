#include <boost/numeric/odeint.hpp>
#include <array>
#include <cstdio>
#include <cmath>
#include <cstring>

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

using ArrayState = std::array<double, 12>;

struct StructState{
	Vec2D  pos [3];
	Vec2D  vel [3];
};

static_assert(sizeof(ArrayState) == sizeof(StructState));

struct System{
	using state = ArrayState;
	double mass[3];
	// ArrayState state;

	// http://adsabs.harvard.edu/full/1967AJ.....72..876S
	void init(ArrayState &ast){
		StructState sst;
		mass[0] = 3.0;
		mass[1] = 4.0;
		mass[2] = 5.0;

		sst.pos [0] = { 1.0,  3.0};
		sst.pos [1] = {-2.0, -1.0};
		sst.pos [2] = { 1.0, -1.0};
		sst.vel [0] = { 0.0,  0.0};
		sst.vel [1] = { 0.0,  0.0};
		sst.vel [2] = { 0.0,  0.0};
		std::memcpy(&ast, &sst, sizeof(StructState));

		this->eng0 = energy(ast);
	}

	void operator()(const ArrayState &aval, ArrayState &ader, double t){
		StructState sval, sder;
		std::memcpy(&sval, &aval, sizeof(StructState));

		auto grav = [](Vec2D ri, Vec2D rj)->Vec2D {
			Vec2D  rij = rj - ri;
			double r2 = rij * rij;
			double r2inv = 1.0 / r2;
			double r3inv = r2inv * std::sqrt(r2inv);

			return r3inv * rij;
		};

		Vec2D f01 = grav(sval.pos[0], sval.pos[1]);
		Vec2D f12 = grav(sval.pos[1], sval.pos[2]);
		Vec2D f20 = grav(sval.pos[2], sval.pos[0]);
		Vec2D f0 = mass[1] * f01 - mass[2] * f20;
		Vec2D f1 = mass[2] * f12 - mass[0] * f01;
		Vec2D f2 = mass[0] * f20 - mass[1] * f12;

		sder.vel[0] = f0;
		sder.vel[1] = f1;
		sder.vel[2] = f2;

		sder.pos[0] = sval.vel[0];
		sder.pos[1] = sval.vel[1];
		sder.pos[2] = sval.vel[2];

		std::memcpy(&ader, &sder, sizeof(StructState));
	}

	double eng0;
	double energy(const ArrayState &aval) const {
		StructState sval;
		std::memcpy(&sval, &aval, sizeof(StructState));

		auto ke = [=](int i)->double {
			return 0.5 * mass[i] * (sval.vel[i] * sval.vel[i]);
		};
		auto pe = [=](int i, int j)->double {
			return -mass[i] *  mass[j] / (sval.pos[j] - sval.pos[i]).abs();
		};

		double kesum = ke(0) + ke(1) + ke(2);
		double pesum = pe(0, 1) + pe(1, 2) + pe(2, 0);

		return kesum + pesum;
	}
};

struct Observer{
	const System &sys;
	void operator()(const ArrayState &s, const double t){
		FILE *fp = stdout;

		double eng = sys.energy(s);
		double derel = (eng - sys.eng0) / sys.eng0;

		fprintf(fp, "%e  %e %e %e %e %e %e  %e %e %e %e %e %e  %e\n",
				t,
				s[0], s[1], s[2], s[3], s[4], s[5],
				s[6], s[7], s[8], s[9], s[10], s[11],
				derel);
													  
	}
};

int main(){
	System sys;
	ArrayState state;

	sys.init(state);

	auto Stepper = boost::numeric::odeint::make_controlled<
		boost::numeric::odeint::runge_kutta_fehlberg78<ArrayState>>(1.e-12, 1.e-10);

	double tend = 50.0;
	boost::numeric::odeint::integrate_adaptive(
			Stepper, sys, state, 0.0, tend, 0.01, Observer{sys});

	return 0;
}
