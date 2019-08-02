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

using ArrayState = std::array<double, 5>;

struct StructState{
	Vec2D  pos;
	Vec2D  vel;
	double time;
};

static_assert(sizeof(ArrayState) == sizeof(StructState));

struct System{
	using state = ArrayState;
	static long neval;

	void init(ArrayState &ast, const double e = 0.5){
		StructState sst;
		sst.pos = {1.0+e, 0.0};
		sst.vel = {0.0, std::sqrt((1.0-e)/(1.0+e))};
		const double r = sst.pos.abs();
		sst.vel = r * sst.vel;
		sst.time = 0.0;

		std::memcpy(&ast, &sst, sizeof(StructState));

		this->eng0 = energy(ast);
	}

	void operator()(const ArrayState &aval, ArrayState &ader, double t){
		StructState sval, sder;
		std::memcpy(&sval, &aval, sizeof(StructState));

		double r = sval.pos.abs();
		double rinv = 1.0 / r;
		double rdotv = sval.pos * sval.vel;
		
		sder.pos = sval.vel;
		sder.vel = rinv * ((rdotv * rinv) * sval.vel - sval.pos);
		sder.time = r;
		std::memcpy(&ader, &sder, sizeof(StructState));

		neval++;
	}

	double eng0;
	double energy(const ArrayState &aval) const {
		StructState sval;
		std::memcpy(&sval, &aval, sizeof(StructState));

		double r = sval.pos.abs();
		double rinv = 1.0 / r;
		double rinv2 = rinv * rinv;

		return -rinv + 0.5 * rinv2 * (sval.vel * sval.vel);
	}
};

long System::neval = 0;

struct Observer{
	const System &sys;
	void operator()(const ArrayState &s, const double tfic){
		FILE *fp = stdout;

		double eng = sys.energy(s);
		double derel = (eng - sys.eng0) / sys.eng0;
		double time = s[4];

		double dt = tfic - tlast;
		tlast = tfic;

		fprintf(fp, "%e  %e %e  %e %e  %e %e %e\n",
				time,
				s[0], s[1], s[2], s[3],
				derel, tfic, dt);
													  
	}
	static double tlast;
};
double Observer::tlast = 0.0;

int main(int ac, char **av){
	double ecc = 0.9;
	if(ac >= 2){
		ecc = atof(av[1]);
	}
	double eabs = 1.e-10;
	double erel = 1.e-10;
	if(ac >= 4){
		eabs = atof(av[2]);
		erel = atof(av[3]);
	}

	System sys;
	ArrayState state;

	sys.init(state, ecc);

#ifndef BSINT
	auto Stepper = boost::numeric::odeint::make_controlled<
		boost::numeric::odeint::runge_kutta_fehlberg78<ArrayState>>(eabs, erel);
#else
	boost::numeric::odeint::bulirsch_stoer<ArrayState> Stepper(eabs, erel);
#endif

	double tend = 10.0;
	boost::numeric::odeint::integrate_adaptive(
			Stepper, sys, state, 0.0, tend, 0.01, Observer{sys});

	fprintf(stderr, "%ld total evaluations\n", System::neval);

	return 0;
}
