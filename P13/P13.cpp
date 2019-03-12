//
// P13.cpp
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
// VE9QRP started the qrpTracker project to fulfill many of the same goals,
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
// - Changed output of method "ascii" to ISO date format

#include "P13.h"


#define DSPMAX_X  480
#define DSPMAX_Y  320

// Converts latitude (Breitengrad) -90..90° / longitude (Laengengrad) -180..180°
// to x/y-coordinates of a map with maxamimum dimension DSPMAX_X * DSPMAX_Y
void latlon2xy(uint16_t &x, uint16_t &y, double lat, double lon)
{
  x = (uint16_t)((((180.0 + lon) / 360.0) * DSPMAX_X) * 0.5);
  y = (uint16_t)(((( 90.0 - lat) / 180.0) * DSPMAX_Y) * 0.5);
}


double
RADIANS(double deg)
{
    return deg * M_PI / 180.0;
}


double
DEGREES(double rad)
{
    return rad * 180.0 / M_PI;
}


//----------------------------------------------------------------------
//     _              ___  _ _______       _      _____ _           
//  __| |__ _ ______ | _ \/ |__ /   \ __ _| |_ __|_   _(_)_ __  ___ 
// / _| / _` (_-<_-< |  _/| ||_ \ |) / _` |  _/ -_)| | | | '  \/ -_)
// \__|_\__,_/__/__/ |_|  |_|___/___/\__,_|\__\___||_| |_|_|_|_\___|
//                                                       
//----------------------------------------------------------------------

static long
fnday(int y, int m, int d)
{
    if (m < 3) {
	m += 12;
	y --;
    }
	
    return (long)(y * YM) + (long)((m+1)*30.6f) + (long)d - 428L;
}


static void
fndate(int &y, int &m, int &d, long dt)
{
    dt += 428L;
    y = (int)((dt-122.1)/365.25);
    dt -= (long)(y*365.25);
    m = (int)(dt / 30.61);
    dt -= (long)(m*30.6);
    m --;
    if (m > 12) {
	m -= 12;
	y++;
    }
    d = dt;
}


P13DateTime::P13DateTime(int year, int month, int day, int h, int m, int s) 
{
    settime(year, month, day, h, m, s);
}


P13DateTime::P13DateTime(const P13DateTime &dt)
{
    DN = dt.DN;
    TN = dt.TN;
}


P13DateTime::P13DateTime()
{
   DN = 0L;
   TN = 0.0;
}


void
P13DateTime::gettime(int &year, int &month, int &day, int &h, int &m, int &s)
{
    fndate(year, month, day, DN);
    double t = TN;
    t *= 24.0;
    h = (int)t;
    t -= h;
    t *= 60;
    m = (int)t;
    t -= m ;
    t *= 60;
    s = (int)t;
}


void
P13DateTime::settime(int year, int month, int day, int h, int m, int s) 
{
    DN = fnday(year, month, day);
    TN = ((double) h + m / 60.0 + s / 3600.0) / 24.0;
}


void
P13DateTime::ascii(char *buf)
{
    int year, mon, day;
    int h, m, s;
	
    gettime(year, mon, day, h, m, s);
    sprintf(buf, "%4d-%02d-%02d %02d:%02d:%02d", year, mon, day, h, m, s);
}


void
P13DateTime::add(double days)
{
    TN += days;
    DN += (long)TN;
    TN -= (long)TN;
}


void
P13DateTime::roundup(double t)
{
    double inc = t - fmod(TN, t);
	
    TN += inc;
    DN += (long)TN;
    TN -= (long)TN;
}

//----------------------------------------------------------------------
//     _              ___  _ ____ ___  _                            
//  __| |__ _ ______ | _ \/ |__ // _ \| |__ ___ ___ _ ___ _____ _ _ 
// / _| / _` (_-<_-< |  _/| ||_ \ (_) | '_ (_-</ -_) '_\ V / -_) '_|
// \__|_\__,_/__/__/ |_|  |_|___/\___/|_.__/__/\___|_|  \_/\___|_|  
//                                                      
//----------------------------------------------------------------------

P13Observer::P13Observer(const char *nm, double lat, double lng, double hgt)
{
    this->name = nm;
    LA = RADIANS(lat);
    LO = RADIANS(lng);
    HT = hgt / 1000.0;

    U[0] = cos(LA) * cos(LO);
    U[1] = cos(LA) * sin(LO);
    U[2] = sin(LA);

    E[0] = -sin(LO);
    E[1] =  cos(LO);
    E[2] =  0.0;

    N[0] = -sin(LA) * cos(LO);
    N[1] = -sin(LA) * sin(LO);
    N[2] =  cos(LA);

    double RP = RE * (1.0 - FL);
    double XX = RE * RE;
    double ZZ = RP * RP;
    double D = sqrt(XX * cos(LA) * cos(LA) + ZZ * sin(LA) * sin(LA));
    double Rx = XX / D + HT;
    double Rz = ZZ / D + HT;

    O[0] = Rx * U[0];
    O[1] = Rx * U[1];
    O[2] = Rz * U[2];

    V[0] = -O[1] * W0;
    V[1] =  O[0] * W0;
    V[2] =  0;
}


//----------------------------------------------------------------------
//     _              ___  _ _______       _       _ _ _ _       
//  __| |__ _ ______ | _ \/ |__ / __| __ _| |_ ___| | (_) |_ ___ 
// / _| / _` (_-<_-< |  _/| ||_ \__ \/ _` |  _/ -_) | | |  _/ -_)
// \__|_\__,_/__/__/ |_|  |_|___/___/\__,_|\__\___|_|_|_|\__\___|
//
//----------------------------------------------------------------------

static double
getdouble(const char *c, int i0, int i1)
{
    char buf[20];
    int i;
	
    for (i=0; i0+i<i1; i++) {
		buf[i] = c[i0+i];
	}
    
	buf[i] = '\0';
    
	return strtod(buf, NULL);
}

static long
getlong(const char *c, int i0, int i1)
{
    char buf[20];
    int i;
	
    for (i=0; i0+i<i1; i++) {
		buf[i] = c[i0+i];
	}
	
    buf[i] = '\0';
	
    return atol(buf);
}

P13Satellite::P13Satellite(const char *nm, const char *l1, const char *l2)
{
    tle(nm, l1, l2);
}

P13Satellite::~P13Satellite()
{
}

void
P13Satellite::tle(const char *nm, const char *l1, const char *l2)
{
    name = nm;

    // direct quantities from the orbital elements

    N  = getlong(l2,  2,  7);
    YE = getlong(l1, 18, 20);
    
	if (YE < 58)
		YE += 2000;
    else
		YE += 1900;

    TE = getdouble(l1, 20, 32);
    M2 = RADIANS(getdouble(l1, 33, 43));

    IN = RADIANS(getdouble(l2, 8, 16));
    RA = RADIANS(getdouble(l2, 17, 25));
    EC = getdouble(l2, 26, 33)/1e7f;
    WP = RADIANS(getdouble(l2, 34, 42));
    MA = RADIANS(getdouble(l2, 43, 51));
    MM = 2.0f * M_PI * getdouble(l2, 52, 63);
    RV = getlong(l2, 63, 68);

    // derived quantities from the orbital elements 

    // convert TE to DE and TE 
    DE = fnday(YE, 1, 0) + (long)TE;
    TE -= (long)TE;
    N0 = MM/86400;
    A_0 = pow(GM/(N0*N0), 1.0/3.0);
    B_0 = A_0*sqrt(1.0-EC*EC);
    PC = RE*A_0/(B_0*B_0);
    PC = 1.5f*J2*PC*PC*MM;
    double CI = cos(IN);
    QD = -PC*CI;
    WD =  PC*(5*CI*CI-1)/2;
    DC = -2*M2/(3*MM);
}


void
P13Satellite::predict(const P13DateTime &dt)
{
    long DN = dt.DN;
    double TN = dt.TN;

    float TEG = DE - fnday(YG, 1, 0) + TE;

    float GHAE = radians(G0) + TEG * WE;
    float MRSE = radians(G0) + TEG * WW + M_PI;
    float MASE = radians(MAS0 + TEG * MASD);

    double T = (double) (DN - DE) + (TN-TE);
    double DT = DC * T / 2.0;
    double KD = 1.0 + 4.0 * DT;
    double KDP = 1.0 - 7.0 * DT;
  
    double M = MA + MM * T * (1.0 - 3.0 * DT);
    double DR = (long) (M / (2.0 * M_PI));
    M -= DR * 2.0 * M_PI;
    double RN = RV + DR;
    double EA = M;

    double DNOM, C_EA, S_EA;

    for (;;) {
		C_EA = cos(EA);
		S_EA = sin(EA);
		DNOM = 1.0 - EC * C_EA;
		double D = (EA - EC * S_EA - M) / DNOM;
		EA -= D ;
		if (fabs(D) < 1e-5)
			break ;
    }

    double A = A_0 * KD;
    double B = B_0 * KD;
    RS = A * DNOM;

    double Vx, Vy;
    double Sx, Sy;
    Sx = A * (C_EA - EC);
    Sy = B * S_EA;

    Vx = -A * S_EA / DNOM * N0;
    Vy =  B * C_EA / DNOM * N0;

    double AP = WP + WD * T * KDP;
    double CW = cos(AP);
    double SW = sin(AP);

    double RAAN = RA + QD * T * KDP;
 
    double CQ = cos(RAAN);
    double SQ = sin(RAAN);

    double CI = cos(IN);
    double SI = sin(IN);

    // CX, CY, and CZ form a 3x3 matrix
    // that converts between orbit coordinates,
    // and celestial coordinates.

    Vec3 CX, CY, CZ;
   
    CX[0] =  CW * CQ - SW * CI * SQ;
    CX[1] = -SW * CQ - CW * CI * SQ;
    CX[2] =  SI * SQ;

    CY[0] =  CW * SQ + SW * CI * CQ;
    CY[1] = -SW * SQ + CW * CI * CQ;
    CY[2] = -SI * CQ;

    CZ[0] = SW * SI;
    CZ[1] = CW * SI;
    CZ[2] = CI;

    // satellite in celestial coords

    SAT[0] = Sx * CX[0] + Sy * CX[1];
    SAT[1] = Sx * CY[0] + Sy * CY[1];
    SAT[2] = Sx * CZ[0] + Sy * CZ[1];

    VEL[0] = Vx * CX[0] + Vy * CX[1];
    VEL[1] = Vx * CY[0] + Vy * CY[1];
    VEL[2] = Vx * CZ[0] + Vy * CZ[1];

    // and in geocentric coordinates

    double GHAA = (GHAE + WE * T);
    double CG = cos(-GHAA);
    double SG = sin(-GHAA);

    S[0] = SAT[0] * CG - SAT[1] * SG;
    S[1] = SAT[0] * SG + SAT[1] * CG;
    S[2] = SAT[2];

    V[0] = VEL[0] * CG - VEL[1]* SG;
    V[1] = VEL[0] * SG + VEL[1]* CG;
    V[2] = VEL[2];
}


void
P13Satellite::LL(double &lat, double &lng)
{
    lat = DEGREES(asin(S[2] / RS));
    lng = DEGREES(atan2(S[1], S[0]));
}


void
P13Satellite::altaz(const P13Observer &obs, double &alt, double &az)
{
    Vec3 R;
	
    R[0] = S[0] - obs.O[0];
    R[1] = S[1] - obs.O[1];
    R[2] = S[2] - obs.O[2];
    double r = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
    R[0] /= r;
    R[1] /= r;
    R[2] /= r;

    double u = R[0] * obs.U[0] + R[1] * obs.U[1] + R[2] * obs.U[2];
    double e = R[0] * obs.E[0] + R[1] * obs.E[1] + R[2] * obs.E[2];
    double n = R[0] * obs.N[0] + R[1] * obs.N[1] + R[2] * obs.N[2];

    az = DEGREES(atan2(e, n));
    
	if (az < 0.0)
		az += 360.0;
    
	alt = DEGREES(asin(u));
}


//----------------------------------------------------------------------
//     _              ___  _ _______           
//  __| |__ _ ______ | _ \/ |__ / __|_  _ _ _  
// / _| / _` (_-<_-< |  _/| ||_ \__ \ || | ' \ 
// \__|_\__,_/__/__/ |_|  |_|___/___/\_,_|_||_|
//                                             
//----------------------------------------------------------------------

Sun::Sun()
{
}

void
P13Sun::predict(const P13DateTime &dt)
{
    long   DN = dt.DN;
    double TN = dt.TN;

    double T = (double) (DN - fnday(YG, 1, 0)) + TN;
    double GHAE = RADIANS(G0) + T * WE;
    double MRSE = RADIANS(G0) + T * WW + M_PI;
    double MASE = RADIANS(MAS0 + T * MASD);
    double TAS = MRSE + EQC1 * sin(MASE) + EQC2 * sin(2.0 * MASE);
    double C, S;

    C = cos(TAS);
    S = sin(TAS);
    SUN[0]=C;
    SUN[1]=S*CNS;
    SUN[2]=S*SNS;
    C = cos(-GHAE); 
    S = sin(-GHAE); 
    H[0]=SUN[0]*C - SUN[1]*S;
    H[1]=SUN[0]*S + SUN[1]*C;
    H[2]=SUN[2];
}

void
P13Sun::LL(double &lat, double &lng)
{
    lat = DEGREES(asin(H[2]));
    lng = DEGREES(atan2(H[1], H[0]));
}

// oops.... this is broken.

void
P13Sun::altaz(const P13Observer &obs, double &alt, double &az)
{
    Vec3 R;

    R[0] = H[0] - obs.O[0];
    R[1] = H[1] - obs.O[1];
    R[2] = H[2] - obs.O[2];
    double r = sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
    R[0] /= r;
    R[1] /= r;
    R[2] /= r;

    double u = R[0] * obs.U[0] + R[1] * obs.U[1] + R[2] * obs.U[2];
    double e = R[0] * obs.E[0] + R[1] * obs.E[1] + R[2] * obs.E[2];
    double n = R[0] * obs.N[0] + R[1] * obs.N[1] + R[2] * obs.N[2];

    az = DEGREES(atan2(e, n));
    
	if (az < 0.0)
		az += 360.0;
    
	alt = DEGREES(asin(u));
}
