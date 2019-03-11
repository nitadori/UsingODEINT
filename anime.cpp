#if 1
#  include <boost/numeric/odeint.hpp>
#else
#  include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#  include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>
#  include <boost/numeric/odeint/stepper/generation/make_controlled.hpp>
#  include <boost/numeric/odeint/stepper/generation/generation_runge_kutta_fehlberg78.hpp>
#endif
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

using ArrayState = std::array<double, 13>;

struct StructState{
	Vec2D  pos [3];
	Vec2D  vel [3];
	double time;
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
		sst.time    = 0.0;
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

		auto Uij = [=](int i, int j)->double {
			auto dr = sval.pos[j] - sval.pos[i];
			auto rinv2 = 1.0 / (dr * dr);;
			// return mass[i] *  mass[j] / (sval.pos[j] - sval.pos[i]).abs();
			return mass[i] *  mass[j] * std::sqrt(rinv2);
		};
		auto Upij = [=](int i, int j)->double {
			auto dr = sval.pos[j] - sval.pos[i];
			auto dv = sval.vel[j] - sval.vel[i];

			auto rdotr = dr * dr;
			auto rdotv = dr * dv;
			auto rinv2 = 1.0 / rdotr;
			auto rinv3 = rinv2 * std::sqrt(rinv2);
			return -(mass[i] *  mass[j]) * rdotv * rinv3;
		};

		Vec2D f01 = grav(sval.pos[0], sval.pos[1]);
		Vec2D f12 = grav(sval.pos[1], sval.pos[2]);
		Vec2D f20 = grav(sval.pos[2], sval.pos[0]);
		Vec2D f0 = mass[1] * f01 - mass[2] * f20;
		Vec2D f1 = mass[2] * f12 - mass[0] * f01;
		Vec2D f2 = mass[0] * f20 - mass[1] * f12;
		
		double U = Uij(0,1) + Uij(1,2) + Uij(2,0);
		double Uprime = Upij(0,1) + Upij(1,2) + Upij(2,0);

		double Uinv = 1.0 / U;
		double Uinv2 = Uinv * Uinv;

		sder.vel[0] = Uinv2 * (f0 - U*Uprime*sval.vel[0]);
		sder.vel[1] = Uinv2 * (f1 - U*Uprime*sval.vel[1]);
		sder.vel[2] = Uinv2 * (f2 - U*Uprime*sval.vel[2]);

		sder.pos[0] = sval.vel[0];
		sder.pos[1] = sval.vel[1];
		sder.pos[2] = sval.vel[2];

		sder.time = Uinv;

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

		return (pesum*pesum)*kesum + pesum;
	}
};

struct Observer{
	const System &sys;

	FILE *gp;
	static constexpr double tick = 0.1;
	double tnext;
	Observer(const System &_sys) : sys(_sys), tnext(0.0)
	{
		gp = popen("gnuplot", "w");
		fprintf(gp, "set terminal gif animate optimize delay 5 size 400,400\n");
		fprintf(gp, "set output 'animation.gif'\n");

		fprintf(gp, "set nokey\n");
		fprintf(gp, "set grid\n");
		fprintf(gp, "set size square\n");
		fprintf(gp, "set xr [-6:6]\n");
		fprintf(gp, "set yr [-6:6]\n");
	}
	~Observer(){
		pclose(gp);
	}

	void operator()(const ArrayState &s, const double tfic){
		FILE *fp = stdout;

		double eng = sys.energy(s);
		double derel = (eng - sys.eng0) / sys.eng0;
		double time = s[12];
		int era = time / 10;

		fprintf(fp, "%e  %e %e %e %e %e %e  %e %e %e %e %e %e  %e %e era%d\n",
				s[12],
				s[0], s[1], s[2], s[3], s[4], s[5],
				s[6], s[7], s[8], s[9], s[10], s[11],
				derel, tfic, era);
													  
		if(time >= tnext){
			fprintf(gp, "set title \"t = %f\"\n", time);
			fprintf(gp, "plot '-' pt 7 ps 1, '-' pt 7 ps 1, '-' pt 7 ps 1, '-' w l lt 1 , '-' w l lt 2, '-' w l lt 3, \n");
			fprintf(gp, "%e, %e\n", s[0], s[1]);
			fprintf(gp, "e\n");
			fprintf(gp, "%e, %e\n", s[2], s[3]);
			fprintf(gp, "e\n");
			fprintf(gp, "%e, %e\n", s[4], s[5]);
			fprintf(gp, "e\n");

			fprintf(gp, "%e %e\n%e %e\ne\n", s[2], s[3], s[4], s[5]);
			fprintf(gp, "%e %e\n%e %e\ne\n", s[4], s[5], s[0], s[1]);
			fprintf(gp, "%e %e\n%e %e\ne\n", s[0], s[1], s[2], s[3]);

			tnext += tick;
		}
	}

};

int main(int ac, char **av){
	double eabs = 1.e-14;
	double erel = 1.e-14;
	if(ac >= 3){
		eabs = atof(av[1]);
		erel = atof(av[2]);
	}

	System sys;
	ArrayState state;

	Observer obs(sys);

	sys.init(state);

	auto Stepper = boost::numeric::odeint::make_controlled<
		boost::numeric::odeint::runge_kutta_fehlberg78<ArrayState>>(eabs, erel);

	double tend = 1950.0;
	boost::numeric::odeint::integrate_adaptive(
			Stepper, sys, state, 0.0, tend, 0.01, std::ref(obs));

	return 0;
}
