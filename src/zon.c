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


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <argp.h>
#include <string.h>
#include <stdlib.h>


const char *argp_program_version =
    PACKAGE_STRING ; 
// e.g. from configure PACKAGE_STRING='zon 2020.11'
const char *argp_program_bug_address =
    PACKAGE_BUGREPORT ; 
// e.g. from configure PACKAGE_BUGREPORT='https://github.com/Aygath/zon'

/* Program documentation. */
static char doc[] =
  "zon [Options...] -- scriptable time output about sun rise and set in UTC";

/* A description of the arguments we accept. */
static char args_doc[] = "-r gives sun rise in iso-timeformat, -s gives sun set, other options are available";

/* The options we understand. */
static struct argp_option options[] = {
  {0,0,0,0, "Options to select what to display" },
  {"rise",     'r', 0,      0,  "Default: Produce the next start/rise time relative to current or provided time" },
  {"set",      's', 0,      0,  "Produce the next end/set time relative to current or provided time" },
  {0,0,0,0, "" },
  {"mid",      'm', 0,      0,  "Produce time exactly between the next rise and set times, i.e. deep midnight or high noon" },
  {"current",  'c', 0,      0,  "Whether sun is up \"+\" or down \"-\"" },
  {0,0,0,0, "Options to specify when and where on earth" },
  {"location", 'l', "+DDMM+DDDMM|+DDMMSS+DDDMMSS", 0,  "Calculate for latitude (+N/-S) and longitude (+E-W) in Degrees, Minutes and Seconds" },
  {"date",     'd', "YYYY-MM-DDTHH:MM+ZZ:zz", 0,  "Calculate for specified iso-formatted time. Defaults to current system UTC date." },
  {0,0,0,0, "Options to format the output. Defaults to iso-format" },
  {"at",       '@', 0,      0,  "Format output as date usable by the at command (in UTC), HH:MM YYYY-MM-DD" },
  {"systemd",  'y', 0,      0,  "Format output as required for systemd-run,  YYYY-MM-DD HH:MM UTC" },
  {"format",   'f', "%H %M %m etc",      0,  "Format output yourself with %H:%M %Y-%m%d %Z etcetera, see strftime() documentation" },
  {0,0,0,0, "" },
  {"verbose",  'v', 0,      0,  "Produce a label or if repeated give all base and calculated data, including date and location" },
  {0,0,0,0, "Options to select the kind of twighlight" },
  {"sun",       0,  0,      0,  "Default: Produce start, ending and duration of visibility of top of sun above horizon, i.e. sunrise and sunset" },
  {0,0,0,0, "" },
  {"civil",     1,  0,      0,  "Produce data about civil twighlight" },
  {0,0,0,0, "" },
  {"nautical",  2,  0,      0,  "Produce data about nautical twighlight" },
  {0,0,0,0, "" },
  {"astronomical",3,0,      0,  "Produce data about nautical twighlight" },
  {0,0,0,0, "" },
  {"angle"     ,5,  "degrees",      0,  "Specify rise/set angle of the centre of the sun" },
  {"rim"       ,4,  0,      0,  "Specify to compensate angle for upper rim of the sun" },
  { 0 }
};
 
/* Used by main to communicate with parse_opt. */
struct arguments
{
  double angle;
  int rise, set, mid, current, verbose, rim;
  char *date, *location, *dateformat;
};

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 0:
	// sun rise (rim of sun above horizon)
      arguments->angle = -35/60;
      arguments->rim =1; 
      break;
    case 1:
      // civil twighlight (center of sun 6 degrees below)
      arguments->angle = -6;
      arguments->rim =0; 
      break;
    case 2:
      // nautical twighlight (center of sun 12 degrees below)
      arguments->angle = -12;
      arguments->rim =0; 
      break;
    case 3:
      // astronimical twighlight (center of sun 18 degrees below)
      arguments->angle = -18;
      arguments->rim =0; 
      break;
    case 4:
      arguments->rim =1; 
      break;
    case 5:
      sscanf(arg,"%lf",&arguments->angle);  
      break;
    case 'c':
      arguments->current = 1;
      break;
    case 'r':
      arguments->rise = 1;
      break;
    case 's':
      arguments->set = 1;
      break;
    case 'm':
      arguments->mid = 1;
      break;
    case '@':
      arguments->dateformat = "%H:%M %Y-%m-%d";
      break;
    case 'I':
      arguments->dateformat = "%Y-%m-%dT%H:%M+00:00";
      break;
    case 'y':
      arguments->dateformat = "%Y-%m-%d %H:%M UTC";
      break;
    case 'f':
      arguments->dateformat = arg;
      break;
    case 'v':
      arguments->verbose += 1;
      break;
    case 'd':
      if ((strlen(arg)!=22) && (strlen(arg)!=21)) argp_usage (state);
      arguments->date = arg;
      break;
    case 'l':
      if (strlen(arg)!=15 && strlen(arg)!=11) argp_usage (state);
      arguments->location = arg;
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

// #static struct argp argp = { 0, 0, 0, doc };
/* Our argp parser. */
// #static struct argp argp = { options, parse_opt, args_doc, doc };
static struct argp argp = { options, parse_opt, 0, doc };

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

int __sunnoonarct__( double d, double lon, double lat,
                  double altit, int upper_limb, double *noon, double *arc );

void sunpos( double d, double *lon, double *r );
void sun_RA_dec( double d, double *RA, double *dec, double *r );
double revolution( double x );
double rev180( double x );
double GMST0( double d );



int main(int argc, char **argv)
{
      struct arguments arguments;
      arguments.current=0;
      arguments.rise=0;
      arguments.set=0;
      arguments.mid=0;
      arguments.verbose=0;
      arguments.angle=-35/60;
      arguments.rim=1;
      arguments.date=NULL;
      arguments.location=NULL;
      arguments.dateformat="%Y-%m-%dT%H:%M+00:00";
      argp_parse (&argp, argc, argv, 0, 0, &arguments);
      if ( ! (arguments.current || arguments.rise || arguments.set || arguments.mid ) ) arguments.rise=1;

      double lon=0,lat=0;
      double tnoon,tarc;
      int    rs;


// Location on earth is our most important parameter for the calculation:
// Successively try: CLI-argument, /etc/zon.conf, TZ environment or /etc/timezone (to inspect /usr/share/zoneinfo.zone1970.tab) or set to 0 Norht, 0 West
// char * getenv (const char *name)
      FILE *fp;
      char *locationstr; locationstr=malloc(sizeof(char)*31);
      char *filename; filename=malloc(sizeof(char)*81);
      if (arguments.location==NULL) { 
	      lon=0; lat=0;
	      strncpy(filename,getenv("HOME"),80);
	      strncat(filename,"/.config/zon.conf",80);
	      fp=fopen(filename,"r");
	      if  (! fp) { 
	        strncpy(filename,"/etc/zon.conf",80);
	        fp=fopen(filename,"r");
	      };
	      if (fp) {
		     if (fgets(locationstr,30,fp)!=0) {
			     locationstr[strlen(locationstr)-1]='\0';
			     // skipt "location=" part:
			     arguments.location=&locationstr[9];
			     if (arguments.verbose >=2 ) printf("Location str from file %s (should be +DDMM[SS]+DDDMM[SS]): %s  len %lu \n",filename,arguments.location,strlen(arguments.location)); 
	             fclose(fp);
		     }
	      }
      } 
      if (arguments.location!=NULL) { 
	 if (arguments.verbose >= 2) printf("Location str (should be +DDMM[SS]+DDDMM[SS]): %s  len %lu \n",arguments.location,strlen(arguments.location)); 
         if ( strlen(arguments.location)==11 ) {
           sscanf(arguments.location, "%5lf%6lf", &lat, &lon); 
           lat = ((int) lat / 100) + (double)((int) lat % 100)/60 ;
           lon = ((int) lon / 100) + (double)((int) lon % 100)/60 ; 
         }
         if ( strlen(arguments.location)==15 ) {
           sscanf(arguments.location, "%7lf%8lf", &lat, &lon); 
           lat = ((int)lat / 10000) + (double)((int)lat / 100 % 100)/60 + (double)((int)lat % 100)/3600 ;
           lon = ((int)lon / 10000) + (double)((int)lon / 100 % 100)/60 + (double)((int)lon % 100)/3600 ; 
         }
      }; 
      if ( arguments.verbose >= 2  ) printf( "Latitude (+ is north) and Longitude (+ is east) decimal values : %+lf %+lf\n", lat, lon );

// moment in time is our next most important reference for the output: 
//      tzset();
      time_t tnow,tbase,trise,tset;
      struct tm base,tmrise,tmset;
      char *datestr; datestr=malloc(sizeof(char)*81);
      double base_d,degree;
      int skipped_days;

      time(&tnow);
      tbase=tnow;
 
// modify base according tot cli-arguments
      gmtime_r(&tbase,&base); 
      int zhours,zmin;
      if ( arguments.date!=NULL) { 
	 if (strlen(arguments.date)==22) 
	   sscanf(arguments.date, "%4d-%2d-%2dT%2d:%2d%3d:%2d", &base.tm_year, &base.tm_mon, &base.tm_mday, &base.tm_hour, &base.tm_min, &zhours, &zmin);  
	 else
	   sscanf(arguments.date, "%4d-%2d-%2dT%2d:%2d%3d%2d", &base.tm_year, &base.tm_mon, &base.tm_mday, &base.tm_hour, &base.tm_min, &zhours, &zmin);  
	 base.tm_year -= 1900; base.tm_mon -= 1; 
	 base.tm_sec = -(zhours*3600 + zmin*60);
	 base.tm_isdst = 0 ;
         tbase = timegm(&base); 
         gmtime_r(&tbase,&base); 
      }; 

      if ( arguments.verbose  >= 2 ) {
	      printf( "system time using local time zone   %s", ctime(&tnow));
              printf( "base time in local time zone        %s", ctime(&tbase));
              printf( "Local Timezone secs                 %ld\n", timezone);
      }; 

// make base data parameters available for calculations:
      /* Compute d of 12h local mean solar time */
      base_d = days_since_2000_Jan_0(base.tm_year+1900,base.tm_mon+1,base.tm_mday) + 0.5 - lon/360.0;

//		int __sunnoonarct_( int year, int month, int day, double lon, double lat,
//                  double altit, int upper_limb, double *noon, double *arc )
	   skipped_days = 0 ; 
           do {
	     rs = __sunnoonarct__(base_d + skipped_days,lon,lat,  arguments.angle, arguments.rim, &tnoon, &tarc ); 
	     trise=tbase + skipped_days*24*60*60;;
	     tmrise= *gmtime(&trise);
	     tmrise.tm_hour = (int)((tnoon-tarc)*60) / 60;
	     tmrise.tm_min  = (int)((tnoon-tarc)*60) % 60;
	     tmrise.tm_sec = 0 ; 
	     tmrise.tm_isdst = 0 ; 
	     trise=timegm(&tmrise); 
	     if (arguments.verbose >= 2 ) if (skipped_days==0 && rs==1) printf("up entire solar day\n"); 
	     ++skipped_days;
           } while  ( !  ( skipped_days>365 || rs==0 && trise > tbase ) ) ; 
	   
	   if (arguments.rise) {
		   gmtime_r(&trise,&tmrise) ;
		   strftime(datestr,80,arguments.dateformat, &tmrise);
		   printf("%s%s\n",datestr,(arguments.verbose>=1)?" rise":"");
	   }
	   
	   skipped_days = 0 ; 
	   do {
	     rs = __sunnoonarct__(base_d + skipped_days,lon,lat, arguments.angle, arguments.rim, &tnoon, &tarc ); 
	     tset=tbase + skipped_days*24*60*60;;
	     tmset= *gmtime(&tset);
	     tmset.tm_hour = (int)((tnoon+tarc)*60) / 60;
	     tmset.tm_min  = (int)((tnoon+tarc)*60) % 60;
	     tmset.tm_sec = 0 ;
	     tmset.tm_isdst = 0 ;
	     tset=timegm(&tmset); 
	     if (arguments.verbose >=2 ) if (skipped_days==0 && rs==-1) printf("down entire solar day\n"); 
	     ++skipped_days;
           } while  ( ! (  skipped_days>365 || rs==0 && tset > tbase ) ) ; 

	   if (arguments.set) {
		   gmtime_r(&tset,&tmset);
		   strftime(datestr,80,arguments.dateformat, &tmset); 
		   printf("%s%s\n",datestr,(arguments.verbose>=1)?" set":"");
	   }

	   if (arguments.current) if (tset<trise)
		   printf("+%s\n",(arguments.verbose >=1)?" up now":"");  
	   else 
		   printf("-%s\n",(arguments.verbose >=1)?" down now":"");  
	
	   if (arguments.mid) {
		   tset = (trise+tset)/2;
		   gmtime_r(&tset,&tmset);
		   strftime(datestr,80,arguments.dateformat, &tmset); 
		   printf("%s%s\n",datestr,(arguments.verbose>=1)?" mid":"");
	   }

      return 0;
}


/* The "workhorse" function for sun south + diurnal arc times */

int __sunnoonarct__( double d, double lon, double lat,
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
{
      double  dd,  /* Days since 2000 Jan 0.0 (negative before) */
      sr,         /* Solar distance, astronomical units */
      sRA,        /* Sun's Right Ascension */
      sdec,       /* Sun's declination */
      sradius,    /* Sun's apparent radius */
//      t,          /* Diurnal arc */
//      tsouth,     /* Time when Sun is at south */
      sidtime;    /* Local sidereal time */

      int rc = 0; /* Return cde from function - usually 0 */

      /* Compute d of 12h local mean solar time */
//      d = days_since_2000_Jan_0(year,month,day) + 0.5 - lon/360.0;

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
