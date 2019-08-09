#include <array>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cassert>

template<
	class    SYS,
	typename T,
	size_t   N>
void check_convergence(
		  SYS &sys, 
		  const T &t, 
		  const T &h, 
		  std::array<T,N> &y,
		  const std::array<T,N> &yref)
{
	using State = std::array<T,N>;
	State f, fmid;
	sys(y, f, t);
	
	State ymid, y1;
	for(size_t i=0; i<N; i++){
		ymid[i] = y[i] + 0.5*h * f[i];
	}

	sys(ymid, fmid, t + 0.5*h);

	for(size_t i=0; i<N; i++){
		y1[i] = y[i] + h * fmid[i];
	}
	 // iteration
	auto break_flag = false;
	for(int n=0; n<100; n++){
		State f1;
		sys(y1, f1, t + h);

		T err_max = 0.0;
		State ynew, ydiff, yerr;
		for(size_t i=0; i<N; i++){
			ynew [i] = y[i] + 0.5*h * (f[i] + f1[i]);
			ydiff[i] = ynew[i] - y1[i];
			y1   [i] = ynew[i];

			err_max = std::max(err_max, std::fabs(ydiff[i]));

			yerr[i] = ynew[i] - yref[i];
		}
		printf("%d %e %e %e %e !conv\n", n, yerr[0], yerr[1], yerr[2], yerr[3]);

		if(break_flag) break;
		if(err_max < 1.e-15) break_flag = true;
	}
}

template<
	class    SYS,
	class    OBS,
	typename T,
	size_t   N>
void integrate(
		  SYS &sys, 
		  OBS obs,
		  const T &t, 
		  const T &h, 
		  std::array<T,N> &y)
{
	using State = std::array<T,N>;

	// 2nd-order initial guess
	State f, fmid;
	sys(y, f, t);
	
	State ymid, y1;
	for(size_t i=0; i<N; i++){
		ymid[i] = y[i] + 0.5*h * f[i];
	}

	sys(ymid, fmid, t + 0.5*h);

	for(size_t i=0; i<N; i++){
		y1[i] = y[i] + h * fmid[i];
	}

	double r = sqrt(y[0]*y[0] + y[1]*y[1]);
	fprintf(stderr, "r = %e\n", r);

	 // iteration
	auto break_flag = false;
	for(int n=0; n<100; n++){
		State f1;
		sys(y1, f1, t + h);

		T err_max = 0.0;
		State ynew, ydiff;
		for(size_t i=0; i<N; i++){
			ynew [i] = y[i] + 0.5*h * (f[i] + f1[i]);
			ydiff[i] = ynew[i] - y1[i];
			y1   [i] = ynew[i];

			err_max = std::max(err_max, std::fabs(ydiff[i]));
		}

#if 0
		fprintf(stderr, "%2d: ", n);
		for(size_t i=0; i<N; i++){
			fprintf(stderr, "%e, ", ydiff[i]);
		}
		fprintf(stderr, "\n");
#endif

		if(true || n > 10){
			fprintf(stderr, "warning n = %d, err_max = %e\n", n, err_max);
		}

		if(break_flag) break;
		if(err_max < 1.e-15) break_flag = true;
	}

	obs(y1, t+h);

	check_convergence(sys, t, h, y, y1);

	// commit
	for(size_t i=0; i<N; i++){
		y[i] = y1[i];
	}

}

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

// #define HARMONIC

struct System{
	using state = ArrayState;
	static long neval;

	void init(ArrayState &ast, const double e = 0.5){
		StructState sst;
#ifndef HARMONIC
		sst.pos = {1.0+e, 0.0};
		sst.vel = {0.0, std::sqrt((1.0-e)/(1.0+e))};
		const double r = sst.pos.abs();
		sst.vel = r * sst.vel;
#else
		sst.pos = {1.0, 0.0};
		sst.vel = {0.0, 1.0};
#endif
		sst.time = 0.0;

		std::memcpy(&ast, &sst, sizeof(StructState));

		this->eng0 = energy(ast);
	}

	void operator()(const ArrayState &aval, ArrayState &ader, double t){
		StructState sval, sder;
		std::memcpy(&sval, &aval, sizeof(StructState));

#ifndef HARMONIC
		double r = sval.pos.abs();
		double rinv = 1.0 / r;
		double rdotv = sval.pos * sval.vel;
		// fprintf(stderr, "rdotv = %e\n", rdotv);
		
		sder.pos = sval.vel;
		sder.vel = rinv * ((rdotv * rinv) * sval.vel - sval.pos);
		sder.time = r;
#else
		sder.pos = sval.vel;
		sder.vel = -1.0 * sval.pos;
		sder.time = 1.0;
#endif
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

		fprintf(fp, "%e %e  %e %e  %e %e  %e %e\n",
				time, tfic,
				s[0], s[1], s[2], s[3],
				derel, dt);
													  
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

	assert(ecc > 0.0);
	assert(ecc <= 1.0);

	System sys;
	ArrayState state;

	sys.init(state, ecc);

	(void)(eabs + erel);

	double h = 0.1;
	double t = 0.0;
	double tend = 3.2;
	while(t < tend){
		integrate(sys, Observer{sys}, t, h, state);
		t += h;
	}


	fprintf(stderr, "%ld total evaluations\n", System::neval);

	return 0;
}
