#include <math.h>

#define rad     57.2957795131e0
#define pi      3.14159265358979e0
#define twopi   6.28318530717959e0

void dcoord(double A0, double B0, double BP,double A1,double B1,double *A2,double *B2)
{
double SB0, CB0, SBP, CBP, SB1,CB1, SB2, CB2;
double SAA, CAA, CBB, SBB, SA2, CA2, TA2O2;

SB0=sin(B0);
CB0=cos(B0);
SBP=sin(BP);
CBP=cos(BP);
SB1=sin(B1);
CB1=cos(B1);
SB2=SBP*SB1 + CBP*CB1*cos(-A1);
CB2=sqrt(1.0e0-pow(SB2,2));
*B2=atan(SB2/CB2);
SAA=sin(-A1)*CB1/CB2;
CAA=(SB1-SB2*SBP)/(CB2*CBP);
CBB=SB0/CBP;
SBB=sin(-A0)*CB0;
SA2=SAA*CBB-CAA*SBB;
CA2=CAA*CBB+SAA*SBB;
TA2O2=0.0e0;															
if(CA2 <= 0.0e0) 
	TA2O2=(1.0e0-CA2)/SA2;  
if(CA2 > 0.0e0) 
	TA2O2=SA2/(1.0e0+CA2);
*A2=2.0e0*atan(TA2O2);
if(*A2 < 0.0e0)
	*A2=(*A2)+6.2831853071795864e0;                                                  
}


// y-Year, m-Month, Day, UTC in hours, RA and Dec of moon
void moon2(unsigned int y, unsigned int m, unsigned int Day, double UT, double lon, double lat, double *RA, double *Dec, double *LST, double *dist, double *Az, double *El)
{
double NN;				//Longitude of ascending node
double i;				//Inclination to the ecliptic
double w;				//Argument of perigee
double a;				//Semi-major axis
double e;				//Eccentricity
double v;				//True anomaly
double d;				//Ephemeris time argument in days
double r;				//Distance to sun, AU
double MM;				//Mean anomaly
double EE;				//Eccentric anomaly
double ecl;				//Obliquity of the ecliptic
double xv, yv;			//x and y coords in ecliptic
double lonecl, latecl;	//Ecliptic long and lat of moon
double xg, yg, zg;		//Ecliptic rectangular coords
double Ms;				//Mean anomaly of sun
double ws;				//Argument of perihelion of sun
double Ls;				//Mean longitude of sun (Ns=0)
double Lm;				//Mean longitude of moon
double DD;				//Mean elongation of moon
double FF;				//Argument of latitude for moon
double xe, ye, ze;		//Equatorial geocentric coords of moon
double mpar;			//Parallax of moon (r_E / d)
double gclat;			//Geocentric latitude
double rho;				//Earth radius factor
double topRA, topDec;	//Topocentric coordinates of Moon
double HA;				//Hour Angle
double GMST0;			//Greenwich Mean Sidereal Time
double A2, B2;
double g; 

//The 'd' in these equations is the decimal day number since 2000 Jan 0.0 UT. This is our time scale. The equation for 'd' follows (y=year, m=month, D=date, UT=UT in hours+decimals):
d = 367 * y - 1.75*(y+(m+9)/12) + 30.5*m +Day - 730530 + UT/24.0e0;     //d = 367*y - 7 * ( y + (m+9)/12 ) / 4 + 275*m/9 + D - 730530 Note: these are ALL integer divisions d = d + UT/24.0
ecl = 23.4393e0 - 3.563e-7 * d;

// Orbital elements for Moon:
NN = 125.1228e0 - 0.0529538083e0 * d;          												
i = 5.1454e0;
w = fmod(318.0634e0 + 0.1643573223e0 * d + 360000.0e0, 360.0e0);          								
a = 60.2666e0;         																	
e = 0.054900e0;        																	
MM = fmod(115.3654e0 + 13.0649929509e0 * d + 360000.0e0, 360.0e0);        								

EE = MM + e*rad*sin(MM/rad) * (1.0e0 + e*cos(MM/rad));         							
EE = EE - (EE - e*rad*sin(EE/rad)-MM) / (1.0e0 - e*cos(EE/rad));
EE = EE - (EE - e*rad*sin(EE/rad)-MM) / (1.0e0 - e*cos(EE/rad));

xv = a * (cos(EE/rad) - e);            														
yv = a * (sqrt(1.0e0-e*e) * sin(EE/rad));

v = fmod(rad*atan2(yv,xv)+720.0e0, 360.0e0);             									
r = sqrt(xv*xv + yv*yv);                														

// Get geocentric position in ecliptic rectangular coordinates:
xg = r * (cos(NN/rad)*cos((v+w)/rad) - sin(NN/rad)*sin((v+w)/rad)*cos(i/rad));    				
yg = r * (sin(NN/rad)*cos((v+w)/rad) + cos(NN/rad)*sin((v+w)/rad)*cos(i/rad));
zg = r * (sin((v+w)/rad)*sin(i/rad)) ;

// Ecliptic longitude and latitude of moon:
lonecl = fmod(rad*atan2(yg/rad,xg/rad)+720.0e0, 360.0e0);
latecl = rad*atan2(zg/rad,sqrt(xg*xg + yg*yg)/rad);

// Now include orbital perturbations:
Ms = fmod(356.0470e0 + 0.9856002585e0 * d + 3600000.0e0, 360.0e0);          								
ws = 282.9404e0 + 4.70935e-5*d;                                 								
Ls = fmod(Ms + ws + 720.0e0, 360.0e0);                          							
Lm = fmod(MM + w + NN+720.0e0, 360.0e0);                        								
DD = fmod(Lm - Ls + 360.0e0, 360.0e0);                         							
FF = fmod(Lm - NN + 360.0e0, 360.0e0);                         						

lonecl = lonecl -1.274e0 * sin((MM-2.0e0*DD)/rad)
                +0.658e0 * sin(2.0e0*DD/rad)
                -0.186e0 * sin(Ms/rad)
                -0.059e0 * sin((2.0e0*MM-2.0e0*DD)/rad)
                -0.057e0 * sin((MM-2.0e0*DD+Ms)/rad)
                +0.053e0 * sin((MM+2.0e0*DD)/rad)
                +0.046e0 * sin((2.0e0*DD-Ms)/rad)
                +0.041e0 * sin((MM-Ms)/rad)
                -0.035e0 * sin(DD/rad)
                -0.031e0 * sin((MM+Ms)/rad)
                -0.015e0 * sin((2.0e0*FF-2.0e0*DD)/rad)
                +0.011e0 * sin((MM-4.0e0*DD)/rad);

latecl = latecl -0.173e0 * sin((FF-2.0e0*DD)/rad)
                -0.055e0 * sin((MM-FF-2.0e0*DD)/rad)
                -0.046e0 * sin((MM+FF-2.0e0*DD)/rad)
                +0.033e0 * sin((FF+2.0e0*DD)/rad)
                +0.017e0 * sin((2.0e0*MM+FF)/rad);

r = 60.36298e0  - 3.27746e0*cos(MM/rad)
                - 0.57994e0*cos((MM-2.0e0*DD)/rad)
                - 0.46357e0*cos(2.0e0*DD/rad)
                - 0.08904e0*cos(2.0e0*MM/rad)
                + 0.03865e0*cos((2.0e0*MM-2.0e0*DD)/rad)
                - 0.03237e0*cos((2.0e0*DD-Ms)/rad)
                - 0.02688e0*cos((MM+2.0e0*DD)/rad)
                - 0.02358e0*cos((MM-2.0e0*DD+Ms)/rad)
                - 0.02030e0*cos((MM-Ms)/rad)
                + 0.01719e0*cos(DD/rad)
                + 0.01671e0*cos((MM+Ms)/rad);

*dist=r*6378.140e0;		//Echo time, seconds

// Geocentric coordinates:
// Rectangular ecliptic coordinates of the moon:
xg = r * cos(lonecl/rad)*cos(latecl/rad);
yg = r * sin(lonecl/rad)*cos(latecl/rad);
zg = r * sin(latecl/rad);

// Rectangular equatorial coordinates of the moon:
xe = xg;                                     													
ye = yg*cos(ecl/rad) - zg*sin(ecl/rad);
ze = yg*sin(ecl/rad) + zg*cos(ecl/rad);

// Right Ascension, Declination:
*RA = fmod(rad*atan2(ye,xe)+360.0e0, 360.0e0);
*Dec = rad*atan2(ze,sqrt(xe*xe + ye*ye));

//Now convert to topocentric system:
mpar=rad*asin(1.0e0/r);
gclat = lat - 0.1924e0*sin(2.0e0*lat/rad);
rho = 0.99883e0 + 0.00167e0*cos(2.0e0*lat/rad);
GMST0 = (Ls + 180.0e0)/15.0e0;
*LST = fmod((GMST0+UT+lon/15.0e0+48.0e0), 24.0e0);		//LST in hours
HA = 15.0e0*(*LST) - (*RA);			//HA in degrees
g = rad*atan(tan(gclat/rad)/cos(HA/rad));
topRA = (*RA) - mpar*rho*cos(gclat/rad)*sin(HA/rad)/cos((*Dec)/rad);
topDec = (*Dec) - mpar*rho*sin(gclat/rad)*sin((g-(*Dec))/rad)/sin(g/rad);

HA = 15.0e0*(*LST) - topRA;		//HA in degrees
if(HA > 180.0e0)
        HA=HA-360.0e0;
if(HA < -180.0e0)
        HA=HA+360.0e0;

dcoord(pi, 0.5e0*pi-lat/rad, lat/rad, HA*twopi/360, topDec/rad, &A2, &B2);

*Az=A2*rad;
*El=B2*rad;
}

void moondoppler(unsigned int y, unsigned int m, unsigned int Day, unsigned int nhr,  unsigned int nmin, unsigned int sec, double nfreq, double lon4, double lat4, double *doppler, double *Az1, double *El1)
{
int i=1;
double vr;			//Radial velocity of moon wrt obs, km/s
double dtopo0;
double rma[10];		//rma -> Vector from Obs to Moon
double rme[10];		//rme -> Vector from Earth center to Moon; 
double rme0[10]; 
double rae[10];		//rae -> Vector from Earth center to Obs; 
double dlat, dlat1, erad1, elev1, dt;
double phi, radps, RA, Dec, LST, dist, Az, El;		// LST -> Locat sidereal time, hours
double UT;		// UT in hours
double f, a, c, arcf, arsf;

UT=(nhr+(nmin/60.0)+(sec/3600.0));

dlat=lat4/rad;		//geocentric
//dlong1=lon4/rad;
elev1=200.0e0;

//-----------------------------geocentric-------------------------------------//
//IAU 1976 flattening f, equatorial radius a
f = 1.0e0/298.257e0;
a = 6378140.0e0;
c = 1.0e0/sqrt(1.0e0 + (-2.0e0 + f)*f*sin(dlat)*sin(dlat));
arcf = (a*c + elev1)*cos(dlat);
arsf = (a*(1.0e0 - f)*(1.0e0 - f)*c + elev1)*sin(dlat);
dlat1 = atan2(arsf,arcf);
erad1 = sqrt(arcf*arcf + arsf*arsf);
erad1 = 0.001e0*(erad1);
//----------------------------------------------------------------------------//

dt=100.0e0;		//For numerical derivative, in seconds

moon2(y,m,Day,(UT-dt/3600.0e0),lon4, lat4, &RA,&Dec,&LST,&dist,&Az,&El);
// Convert to rectangular coords
rme0[1]=dist*cos(Dec/rad)*cos(RA/rad);          												
rme0[2]=dist*cos(Dec/rad)*sin(RA/rad);
rme0[3]=dist*sin(Dec/rad);

moon2(y,m,Day,UT,lon4,lat4,&RA,&Dec,&LST,&dist,&Az,&El);
*Az1=Az;
*El1=El;
rme[1]=dist*cos(Dec/rad)*cos(RA/rad);		//Convert to rectangular coords
rme[2]=dist*cos(Dec/rad)*sin(RA/rad);
rme[3]=dist*sin(Dec/rad);

phi=LST*twopi/24.0e0;         																
// Gencentric numbers used here!
rae[1]=erad1*cos(dlat1)*cos(phi);
rae[2]=erad1*cos(dlat1)*sin(phi);
rae[3]=erad1*sin(dlat1);

radps=twopi/(86400.0e0/1.002737909e0);      																		    

rae[4]= -rae[2]*radps;		//Vel of Obs wrt Earth center
rae[5]= rae[1]*radps;
rae[6]= 0.0e0 ;

for(i=1;i<4;i++)
{rme[i+3]=(rme[i]-rme0[i])/dt;
 rma[i]=rme[i]-rae[i];
 rma[i+3]=rme[i+3]-rae[i+3];}

dtopo0=sqrt(pow(rma[1],2) + pow(rma[2],2) + pow(rma[3],2));		//Get topocentric coords
vr = ( rma[1]*rma[4] + rma[2]*rma[5] + rma[3]*rma[6] ) / dtopo0;		//Radial velocity of moon wrt obs, km/s

*doppler=-nfreq * vr/2.99792458e5 * pow(10,6) * 2;		//doppler=-freq*vr/2.99792458e5 ; freq = freq*10^6 to transform in Hz, *2 doppler in both directions
}

