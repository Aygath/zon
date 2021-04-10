// zon.c
// Last modified march 2021
// This program will print the date/time of the next sunrise and/or sunset in a format that
// is convenient for scripting and scheduling with e.g. the at command and the systemd command.
// Example echo myjob.sh | TZ=UTC at $(zon -@t)
// Example echo Next sunrise in UTC will be at $(zon -r)
// Example echo Next sunrise in in local time will be at $(date -d $(zon -r) ) 
//
// The actual calculations are copied from the programs SUNRISET.C written by 
// Paul Schlyter and released to the public domain in december 1992. See comment below.
//
/* +++Date last modified: 05-Jul-1997 */
/* Updated comments, 05-Aug-2013 */
 
/*

SUNRISET.C - computes Sun rise/set times, start/end of twilight, and
             the length of the day at any date and latitude

Written as DAYLEN.C, 1989-08-16

Modified to SUNRISET.C, 1992-12-01

(c) Paul Schlyter, 1989, 1992

Released to the public domain by Paul Schlyter, December 1992

*/


/* A macro to compute the number of days elapsed since 2000 Jan 0.0 */
/* (which is equal to 1999 Dec 31, 0h UT)                           */

#define days_since_2000_Jan_0(y,m,d) \
    (367L*(y)-((7*((y)+(((m)+9)/12)))/4)+((275*(m))/9)+(d)-730530L)

/* Some conversion factors between radians and degrees */

#ifndef PI
 #define PI        3.1415926535897932384
#endif

#define RADEG     ( 180.0 / PI )
#define DEGRAD    ( PI / 180.0 )

/* The trigonometric functions in degrees */

#define sind(x)  sin((x)*DEGRAD)
#define cosd(x)  cos((x)*DEGRAD)
#define tand(x)  tan((x)*DEGRAD)

#define atand(x)    (RADEG*atan(x))
#define asind(x)    (RADEG*asin(x))
#define acosd(x)    (RADEG*acos(x))
#define atan2d(y,x) (RADEG*atan2(y,x))



/* Function prototypes */

void sunpos( double d, double *lon, double *r );
void sun_RA_dec( double d, double *RA, double *dec, double *r );
double revolution( double x );
double rev180( double x );
double GMST0( double d );
void moonpos( double d, double *lon, double *lat, double *r );
void mon_RA_dec( double d, double *RA, double *dec, double *r );



int __sunnoonarct__( int year, int month, int day, double lon, double lat,
                  double altit, int upper_limb, double *noon, double *arc )
/***************************************************************************/
/* Note: year,month,date = calendar date, 1801-2099 only.             */
/*       Eastern longitude positive, Western longitude negative       */
/*       Northern latitude positive, Southern latitude negative       */
/*       The longitude value IS critical in this function!            */
/*       altit = the altitude which the Sun should cross              */
/*               Set to -35/60 degrees for rise/set, -6 degrees       */
/*               for civil, -12 degrees for nautical and -18          */
/*               degrees for astronomical twilight.                   */
/*         upper_limb: non-zero -> upper limb, zero -> center         */
/*               Set to non-zero (e.g. 1) when computing rise/set     */
/*               times, and to zero when computing start/end of       */
/*               twilight.                                            */
/*        *rise = where to store the rise time                        */
/*        *set  = where to store the set  time                        */
/*                Both times are relative to the specified altitude,  */
/*                and thus this function can be used to compute       */
/*                various twilight times, as well as rise/set times   */
/* Return value:  0 = sun rises/sets this day, times stored at        */
/*                    *trise and *tset.                               */
/*               +1 = sun above the specified "horizon" 24 hours.     */
/*                    *trise set to time when the sun is at south,    */
/*                    minus 12 hours while *tset is set to the south  */
/*                    time plus 12 hours. "Day" length = 24 hours     */
/*               -1 = sun is below the specified "horizon" 24 hours   */
/*                    "Day" length = 0 hours, *trise and *tset are    */
/*                    both set to the time when the sun is at south.  */
/*                                                                    */
/**********************************************************************/
/* This procedure calculates a time on a given date, i.e. it is       */
/* independent of an exact time during that day.                      */
{
      
      double sr,  /* Solar distance, astronomical units */
      sRA,        /* Sun's Right Ascension */
      sdec,       /* Sun's declination */
      sradius,    /* Sun's apparent radius */
      d,          /* Days since 2000 Jan 0.0 (negative before) */
//      t,          /* Diurnal arc */
//      tsouth,     /* Time when Sun is at south */
      sidtime;    /* Local sidereal time */

      int rc = 0; /* Return cde from function - usually 0 */

      /* Compute d of 12h local mean solar time */
      d = days_since_2000_Jan_0(year,month,day) + 0.5 - lon/360.0;

      /* Compute the local sidereal time of this moment */
      sidtime = revolution( GMST0(d) + 180.0 + lon );

      /* Compute Sun's RA, Decl and distance at this moment */
      sun_RA_dec( d, &sRA, &sdec, &sr );

      /* Compute time when Sun is at south - in hours UT */
      *noon = 12.0 - rev180(sidtime - sRA)/15.0;

      /* Compute the Sun's apparent radius in degrees */
      sradius = 0.2666 / sr;

      /* Do correction to upper limb, if necessary */
      if ( upper_limb )
            altit -= sradius;

      /* Compute the diurnal arc that the Sun traverses to reach */
      /* the specified altitude altit: */
      {
            double cost;
            cost = ( sind(altit) - sind(lat) * sind(sdec) ) /
                  ( cosd(lat) * cosd(sdec) );
            if ( cost >= 1.0 )
                  rc = -1, *arc = 0.0;       /* Sun always below altit */
            else if ( cost <= -1.0 )
                  rc = +1, *arc = 12.0;      /* Sun always above altit */
            else
                  *arc = acosd(cost)/15.0;   /* The diurnal arc, hours */
      }
      /* Store rise and set times - in hours UT */
//      *trise = tsouth - t;
//      *tset  = tsouth + t;

      return rc;
}  /* __sunriset__ */



/* The "workhorse" function for sun rise/set times */

int __sunriset__( int year, int month, int day, double lon, double lat,
                  double altit, int upper_limb, double *trise, double *tset )
/***************************************************************************/
/* Note: year,month,date = calendar date, 1801-2099 only.             */
/*       Eastern longitude positive, Western longitude negative       */
/*       Northern latitude positive, Southern latitude negative       */
/*       The longitude value IS critical in this function!            */
/*       altit = the altitude which the Sun should cross              */
/*               Set to -35/60 degrees for rise/set, -6 degrees       */
/*               for civil, -12 degrees for nautical and -18          */
/*               degrees for astronomical twilight.                   */
/*         upper_limb: non-zero -> upper limb, zero -> center         */
/*               Set to non-zero (e.g. 1) when computing rise/set     */
/*               times, and to zero when computing start/end of       */
/*               twilight.                                            */
/*        *rise = where to store the rise time                        */
/*        *set  = where to store the set  time                        */
/*                Both times are relative to the specified altitude,  */
/*                and thus this function can be used to compute       */
/*                various twilight times, as well as rise/set times   */
/* Return value:  0 = sun rises/sets this day, times stored at        */
/*                    *trise and *tset.                               */
/*               +1 = sun above the specified "horizon" 24 hours.     */
/*                    *trise set to time when the sun is at south,    */
/*                    minus 12 hours while *tset is set to the south  */
/*                    time plus 12 hours. "Day" length = 24 hours     */
/*               -1 = sun is below the specified "horizon" 24 hours   */
/*                    "Day" length = 0 hours, *trise and *tset are    */
/*                    both set to the time when the sun is at south.  */
/*                                                                    */
/**********************************************************************/
{
      double  d,  /* Days since 2000 Jan 0.0 (negative before) */
      sr,         /* Solar distance, astronomical units */
      sRA,        /* Sun's Right Ascension */
      sdec,       /* Sun's declination */
      sradius,    /* Sun's apparent radius */
      t,          /* Diurnal arc */
      tsouth,     /* Time when Sun is at south */
      sidtime;    /* Local sidereal time */

      int rc = 0; /* Return cde from function - usually 0 */

      /* Compute d of 12h local mean solar time */
      d = days_since_2000_Jan_0(year,month,day) + 0.5 - lon/360.0;

      /* Compute the local sidereal time of this moment */
      sidtime = revolution( GMST0(d) + 180.0 + lon );

      /* Compute Sun's RA, Decl and distance at this moment */
      sun_RA_dec( d, &sRA, &sdec, &sr );

      /* Compute time when Sun is at south - in hours UT */
      tsouth = 12.0 - rev180(sidtime - sRA)/15.0;

      /* Compute the Sun's apparent radius in degrees */
      sradius = 0.2666 / sr;

      /* Do correction to upper limb, if necessary */
      if ( upper_limb )
            altit -= sradius;

      /* Compute the diurnal arc that the Sun traverses to reach */
      /* the specified altitude altit: */
      {
            double cost;
            cost = ( sind(altit) - sind(lat) * sind(sdec) ) /
                  ( cosd(lat) * cosd(sdec) );
            if ( cost >= 1.0 )
                  rc = -1, t = 0.0;       /* Sun always below altit */
            else if ( cost <= -1.0 )
                  rc = +1, t = 12.0;      /* Sun always above altit */
            else
                  t = acosd(cost)/15.0;   /* The diurnal arc, hours */
      }

      /* Store rise and set times - in hours UT */
      *trise = tsouth - t;
      *tset  = tsouth + t;

      return rc;
}  /* __sunriset__ */




/* This function computes the Sun's position at any instant */

void sunpos( double d, double *lon, double *r )
/******************************************************/
/* Computes the Sun's ecliptic longitude and distance */
/* at an instant given in d, number of days since     */
/* 2000 Jan 0.0.  The Sun's ecliptic latitude is not  */
/* computed, since it's always very near 0.           */
/******************************************************/
{
      double M,         /* Mean anomaly of the Sun */
             w,         /* Mean longitude of perihelion */
                        /* Note: Sun's mean longitude = M + w */
             e,         /* Eccentricity of Earth's orbit */
             E,         /* Eccentric anomaly */
             x, y,      /* x, y coordinates in orbit */
             v;         /* True anomaly */

      /* Compute mean elements */
      M = revolution( 356.0470 + 0.9856002585 * d );
      w = 282.9404 + 4.70935E-5 * d;
      e = 0.016709 - 1.151E-9 * d;

      /* Compute true longitude and radius vector */
      E = M + e * RADEG * sind(M) * ( 1.0 + e * cosd(M) );
            x = cosd(E) - e;
      y = sqrt( 1.0 - e*e ) * sind(E);
      *r = sqrt( x*x + y*y );              /* Solar distance */
      v = atan2d( y, x );                  /* True anomaly */
      *lon = v + w;                        /* True solar longitude */
      if ( *lon >= 360.0 )
            *lon -= 360.0;                   /* Make it 0..360 degrees */
}



void sun_RA_dec( double d, double *RA, double *dec, double *r )
/******************************************************/
/* Computes the Sun's equatorial coordinates RA, Decl */
/* and also its distance, at an instant given in d,   */
/* the number of days since 2000 Jan 0.0.             */
/******************************************************/
{
      double lon, obl_ecl, x, y, z;

      /* Compute Sun's ecliptical coordinates */
      sunpos( d, &lon, r );

      /* Compute ecliptic rectangular coordinates (z=0) */
      x = *r * cosd(lon);
      y = *r * sind(lon);

      /* Compute obliquity of ecliptic (inclination of Earth's axis) */
      obl_ecl = 23.4393 - 3.563E-7 * d;

      /* Convert to equatorial rectangular coordinates - x is unchanged */
      z = y * sind(obl_ecl);
      y = y * cosd(obl_ecl);

      /* Convert to spherical coordinates */
      *RA = atan2d( y, x );
      *dec = atan2d( z, sqrt(x*x + y*y) );

}  /* sun_RA_dec */

void moonpos( double d, double *lon, double *lat, double *r )
/******************************************************/
/* Computes the Moon's ecliptic longitude, latitude   */
/* and distance at an instant given in d, number of   */
/* days since 2000 Jan 0.0.                           */
/* Adapted from sunpos by Aygath/github               */
/******************************************************/
// The Moon's position, as computed, is geocentric, i.e. as seen by an imaginary observer at the center of the Earth
{
      double M,         /* Mean anomaly of the Moon */
             w,         /* Mean longitude of perihelion */
	                /* Sun's mean longitude Ls = Ms + ws */
                        /* Note: Moon's mean longitude Lm = N + M + w */
             e,         /* Eccentricity of Earth's orbit */
	     N,		/* (Long asc. node) */
	     i,         /* Inclination */
	     a,		/* Mean distance */
             E0,E,      /* Eccentric anomaly */
             x, y,      /* x, y coordinates in orbit */
	     xeclip,
	     yeclip,
	     zeclip,    /* ecliptic coordinates */
             v;         /* True anomaly */

//     printf ("d=%f\n",d);
      /* Compute mean elements */
      M = revolution( 115.3654 + 13.0649929509 * d );
      w = revolution(318.0634 + 0.1643573223 * d);
      e = 0.054900;
      N = revolution(125.1228 - 0.0529538083  * d);
      i = 5.1454;
      a = 60.2666;


//      printf("M=%f  w=%f  N=%f\n",M,w,N);
      /* Compute true longitude and radius vector */
      E0 = M + e * RADEG * sind(M) * ( 1.0 + e * cosd(M) );
      E = E0 - (E0 - e * RADEG * sind(E0)-M)/(1.0 + e * cosd(E0));
//	      printf("E=%f\n",E);
      while (abs(E0-E)>0.005) {  
	      E0=E;
              E = E0 - (E0 - ((e * RADEG * sind(E0)-M)/(1.0 + e * cosd(E0))));
//	      printf("E=%f\n",E);
      } ;

      x = a * (cosd(E) - e);
      y = a * sqrt( 1.0 - e*e ) * sind(E);

      *r = sqrt( x*x + y*y );              /* Moon's distance */
      v =  revolution( atan2d( y, x )) ;                  /* True anomaly */
//      v =  259.8605;                  /* True anomaly */
//      printf("x=%f  y=%f   r=%f  v=%f \n",x,y,*r,v);

      xeclip = *r * ( cosd(N) * cosd(v+w) - sind(N) * sind(v+w) * cosd(i) );
      yeclip = *r * ( sind(N) * cosd(v+w) + cosd(N) * sind(v+w) * cosd(i) );
      zeclip = *r * sind(v+w) * sind(i);
//printf("xec=%f   yec=%f   zec=%f\n",xeclip,yeclip,zeclip);

      *lon =revolution atan2d( yeclip, xeclip );                        /* True moon's longitude */
      *lat = atan2d( zeclip, sqrt( xeclip*xeclip + yeclip*yeclip ) );       /* Moon's latitude */
      /* in fact optional : */ 
      *r = sqrt( xeclip*xeclip + yeclip*yeclip + zeclip*zeclip );

      /* now add in perturbations cause by the sun */ 
      /* Compute mean element for the Sun */
      double Ms = revolution( 356.0470 + 0.9856002585 * d ); /* Mean anomaly of the Sun  */
      /* Save some funcamental arguments for convenience */
      double D = revolution( (N + w + M) - (Ms + (282.9404 + 4.70935E-5 * d) )); /* Moon's mean elongation */
      double F = revolution( w + M ); /* Moon's argument of latitude */
// printf("M=%f  Ms=%f  D=%f  F=%f\n",M,Ms,D,F);
      *lon += (-1.274 * sind(M - 2*D))
	      + (0.658 * sind(2*D))
	      + (-0.186 * sind(Ms))
	      + (-0.059 * sind(2*M - 2*D))
	      + (-0.057 * sind(M - 2*D + Ms))
	      + (0.053 * sind(M + 2*D))
	      + (0.046 * sind(2*D - Ms))
	      + (0.041 * sind(M - Ms))
	      + (-0.035 * sind(D))
	      + (-0.031 * sind(M + Ms))
	      + (-0.015 * sind(2*F - 2*D))
	      + (0.011 * sind(M - 4*D))
	      ;
      *lat += (-0.173 * sind(F - 2*D))
	      + (-0.055 * sind(M - F - 2*D))
	      + (-0.046 * sind(M + F - 2*D))
	      + (0.033 * sind(F + 2*D))
	      + (0.017 * sind(2*M + F)) /* why does this term not relate to the sun ? */
	      ;
      *r   += (-0.58 * cosd(M - 2*D))
	      + (-0.46 * cosd(2*D))
	      ;
}


void moon_RA_dec( double d, double *RA, double *dec, double *r )
/******************************************************/
/* Computes the Sun's equatorial coordinates RA, Decl */
/* and also its distance, at an instant given in d,   */
/* the number of days since 2000 Jan 0.0.             */
/******************************************************/
{
      double lon, lat, obl_ecl, x, y, z, y_eq;

      /* Compute Sun's ecliptical coordinates */
      moonpos( d, &lon, &lat, r );

      /* Compute ecliptic rectangular coordinates (z=0) */
      x = cosd(lon)*cosd(lat);
      y = sind(lon)*cosd(lat);
      z = sind(lat);
      /* Compute obliquity of ecliptic (inclination of Earth's axis) */
      obl_ecl = 23.4393 - 3.563E-7 * d;

      /* Convert to equatorial rectangular coordinates - x is unchanged */
      y_eq = y * cosd(obl_ecl) - z * sind(obl_ecl);
      z = y * sind(obl_ecl) + z * cosd(obl_ecl);

      /* Convert to GEOCENTRIC spherical coordinates */
      *RA = revolution( atan2d( y_eq, x ) );
      *dec = atan2d( z, sqrt(x*x + y_eq*y_eq) );

      /* Now make it TOPOCENTRIC, because the moon is "close" to Earth */
      double mpar, gclat, rho,HA, UT, g;

      /* Moon's parallax, i.e. disk of Earth as seen from Moon */
      mpar = asind(1.0 / *r);
      gclat = lat - 0.1924 * sind(2.0*lat);
      rho   = 0.99833 + 0.00167 * cosd(2.0*lat);
//      HA = (GMST0(d)+ UT + lon/15.0 ) - *RA;
      g = atand(tand(gclat)/cosd(HA));
      *RA = *RA - mpar * rho * cosd(gclat) * sind(HA)/cosd(*dec);
      *dec = *dec - mpar * rho * sind(gclat) * sind(g - *dec)/sind(g);

}  /* sun_RA_dec */
/******************************************************************/
/* This function reduces any angle to within the first revolution */
/* by subtracting or adding even multiples of 360.0 until the     */
/* result is >= 0.0 and < 360.0                                   */
/******************************************************************/

#define INV360    ( 1.0 / 360.0 )

double revolution( double x )
/*****************************************/
/* Reduce angle to within 0..360 degrees */
/*****************************************/
{
      return( x - 360.0 * floor( x * INV360 ) );
}  /* revolution */

double rev180( double x )
/*********************************************/
/* Reduce angle to within +180..+180 degrees */
/*********************************************/
{
      return( x - 360.0 * floor( x * INV360 + 0.5 ) );
}  /* revolution */


/*******************************************************************/
/* This function computes GMST0, the Greenwich Mean Sidereal Time  */
/* at 0h UT (i.e. the sidereal time at the Greenwhich meridian at  */
/* 0h UT).  GMST is then the sidereal time at Greenwich at any     */
/* time of the day.  I've generalized GMST0 as well, and define it */
/* as:  GMST0 = GMST - UT  --  this allows GMST0 to be computed at */
/* other times than 0h UT as well.  While this sounds somewhat     */
/* contradictory, it is very practical:  instead of computing      */
/* GMST like:                                                      */
/*                                                                 */
/*  GMST = (GMST0) + UT * (366.2422/365.2422)                      */
/*                                                                 */
/* where (GMST0) is the GMST last time UT was 0 hours, one simply  */
/* computes:                                                       */
/*                                                                 */
/*  GMST = GMST0 + UT                                              */
/*                                                                 */
/* where GMST0 is the GMST "at 0h UT" but at the current moment!   */
/* Defined in this way, GMST0 will increase with about 4 min a     */
/* day.  It also happens that GMST0 (in degrees, 1 hr = 15 degr)   */
/* is equal to the Sun's mean longitude plus/minus 180 degrees!    */
/* (if we neglect aberration, which amounts to 20 seconds of arc   */
/* or 1.33 seconds of time)                                        */
/*                                                                 */
/*******************************************************************/

double GMST0( double d )
{
      double sidtim0;
      /* Sidtime at 0h UT = L (Sun's mean longitude) + 180.0 degr  */
      /* L = M + w, as defined in sunpos().  Since I'm too lazy to */
      /* add these numbers, I'll let the C compiler do it for me.  */
      /* Any decent C compiler will add the constants at compile   */
      /* time, imposing no runtime or code overhead.               */
      sidtim0 = revolution( ( 180.0 + 356.0470 + 282.9404 ) +
                          ( 0.9856002585 + 4.70935E-5 ) * d );
      return sidtim0;
}  /* GMST0 */

void EqAz( double RA, double DEC, struct tm tnow , double lon, double lat, double *azimuth, double *altitude)
/* This converts the RA:decl angle to azimuth:altitude angle                        */
/* You must specify your reference frame for azi:alt by supplying you location plus */
/* the time (d) at which you want to observe RA:DEC */
{
      double HA, x,y,z, d, xhor,yhor,zhor,LST;
      /* Compute d of local mean solar time */
      d = days_since_2000_Jan_0(tnow.tm_year+1900,tnow.tm_mon+1,tnow.tm_mday) + tnow.tm_hour/24.0 + tnow.tm_min/(24*60.0) - lon/360.0;
      /* Compute the local sidereal time of this moment */
      LST = revolution( GMST0(d) + 180.0 + lon );
      HA = LST - RA;
      x = cosd(HA) * cosd(DEC);
      y = sind(HA) * cosd(DEC);
      z = sind(DEC);

      xhor = x * sind(lat) - z * cosd(lat);
      yhor = y;
      zhor = x * cosd(lat) + z * sind(lat);

      *azimuth  = atan2d( yhor, xhor ) + 180;
      *altitude = asind( zhor )       ;
}

