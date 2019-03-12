//
// Plan13.cpp
//
// An implementation of Plan13 in C++ by Mark VandeWettering
//
// Plan13 is an algorithm for satellite orbit prediction first formulated
// by James Miller G3RUH.  I learned about it when I saw it was the basis 
// of the PIC based antenna rotator project designed by G6LVB.
//
// http://www.g6lvb.com/Articles/LVBTracker2/index.htm
//
// I ported the algorithm to Python, and it was my primary means of orbit
// prediction for a couple of years while I operated the "Easy Sats" with 
// a dual band hand held and an Arrow antenna.
//
// I've long wanted to redo the work in C++ so that I could port the code
// to smaller processors including the Atmel AVR chips.  Bruce Robertson,
// VE9QRP started the qrpTracker project to fufill many of the same goals,
// but I thought that the code could be made more compact and more modular,
// and could serve not just the embedded targets but could be of more
// use for more general applications.  And, I like the BSD License a bit
// better too.
//
// So, here it is!
//

// 26.02.2019 dl9sec@gmx.net:
// - Renamed class DateTime to P13DateTime because of potential conflict
//   with RTClib, which uses the same class name.
// - Renamed all classes to P13... for consistency and clarification of
//   affiliation.
// - Some code beautification.
// - Fixed missing semicolon at destructor of P13DateTime.
// - Added comments for constants and TLE data.
// - Updated sidereal and solar data to values valid until ~2030

//----------------------------------------------------------------------

#include "WProgram.h"

// here are a bunch of constants that will be used throughout the 
// code, but which will probably not be helpful outside.

static const double RE = 6378.137f;				// WGS-84 Earth ellipsoid
static const double FL = 1.0f/298.257224f;		// -"-
static const double GM = 3.986E5f;				// Earth's Gravitational constant km^3/s^2
static const double J2 = 1.08263E-3f;			// 2nd Zonal coeff, Earth's Gravity Field
static const double YM = 365.25f;				// Mean Year,     days
static const double YT = 365.2421874f;			// Tropical year, days
static const double WW = 2.0f*M_PI/YT;			// Earth's rotation rate, rads/whole day
static const double WE = 2.0f*M_PI+ WW;			// Earth's rotation rate, radians/day 
static const double W0 = WE/86400.0f;			// Earth's rotation rate, radians/sec

// Sidereal and Solar data. Rarely needs changing. Valid to year ~2015
/*
static const double YG = 2000.0f;				// GHAA, Year YG, Jan 0.0
static const double G0 = 98.9821f;				// -"-
static const double MAS0 = 356.0507f;			// MA Sun and rate, deg, deg/day
static const double MASD = 0.98560028f;			// -"-
static const double EQC1 = 0.03342;				// Sun's Equation of centre terms
static const double EQC2 = 0.00035;				// -"-
static const double INS = radians(23.4393f);	// Sun's inclination
static const double CNS = cos(INS);				// -"-
static const double SNS = sin(INS);				// -"-
*/

// Sidereal and Solar data. Rarely needs changing. Valid to year ~2030
static const double YG = 2014.0f;				// GHAA, Year YG, Jan 0.0
static const double G0 = 99.5828f;				// -"-
static const double MAS0 = 356.4105f;			// MA Sun and rate, deg, deg/day
static const double MASD = 0.98560028f;			// -"-
static const double EQC1 = 0.03340;				// Sun's Equation of centre terms
static const double EQC2 = 0.00035;				// -"-
static const double INS = radians(23.4375f);	// Sun's inclination
static const double CNS = cos(INS);				// -"-
static const double SNS = sin(INS);				// -"-

//----------------------------------------------------------------------

// the original BASIC code used three variables (e.g. Ox, Oy, Oz) to
// represent a vector quantity.  I think that makes for slightly more
// obtuse code, so I going to collapse them into a single variable 
// which is an array of three elements

typedef double Vec3[3];

//----------------------------------------------------------------------

class P13DateTime {
public:
	long   DN;
 	double TN;
	
   	P13DateTime(int year, int month, int day, int h, int m, int s);
	P13DateTime(const P13DateTime &);
	P13DateTime();
	~P13DateTime() { };
	void add(double);
   	void settime(int year, int month, int day, int h, int m, int s);
    void gettime(int& year, int& mon, int& day, int& h, int& m, int& s);
    void ascii(char *);
	void roundup(double);
};


//----------------------------------------------------------------------

class P13Observer {
public:
    const char *name;
    double LA;
    double LO;
    double HT;
    Vec3 U, E, N, O, V;
    
    P13Observer(const char *, double, double, double);
    ~P13Observer() { };
};


//----------------------------------------------------------------------

class P13Satellite { 
  	long   N;
	long   YE;		// Epoch Year    			year
    long   DE;		// Epoch Fraction of day
	double TE;		// Epoch time    			days
	double IN;		// Inclination   			deg
	double RA;		// R.A.A.N.      			deg
	double EC;		// Eccentricity  			 -
	double WP;		// Arg perigee   			deg
	double MA;		// Mean anomaly  			deg
	double MM;		// Mean motion   			rev/d
	double M2;		// Decay Rate    			rev/d/d
	double RV;		// Orbit number  			 -
	double ALON;	// Sat attitude				deg
	double ALAT;	// Sat attitude				deg

	// these values are stored, but could be calculated on the fly
	// during calls to predict() 
	// classic space/time tradeoff

    double N0, A_0, B_0;
    double PC;
    double QD, WD, DC;
    double RS;

public:
    const char *name;
	Vec3 SAT, VEL;		// celestial coordinates
    Vec3 S, V; 			// geocentric coordinates
 
	P13Satellite() { };
	P13Satellite(const char *name, const char *l1, const char *l2);
	~P13Satellite();
    void tle(const char *name, const char *l1, const char *l2);
    void predict(const P13DateTime &dt);
 	void LL(double &lat, double &lng);
	void altaz(const P13Observer &obs, double &alt, double &az);
};


//----------------------------------------------------------------------

class P13Sun {
public:
	Vec3 SUN, H;
	P13Sun();
	~P13Sun() { };
    void predict(const P13DateTime &dt);
 	void LL(double &lat, double &lng);
	void altaz(const P13Observer &obs, double &alt, double &az);
};
