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
#include <regex.h>
#include "astro.c"


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
  {"rise",     'r', 0,      0,  "Produce the next start/rise time relative to current or provided time" },
  {"set",      's', 0,      0,  "Produce the next end/set time relative to current or provided time" },
  {0,0,0,0, "" },
  {"mid",      'm', 0,      0,  "Produce time exactly between the next rise and set times, i.e. deep midnight or high noon" },
  {"current",  'c', 0,      0,  "Whether sun is up \"+\" or down \"-\". Default if no display type selected." },
  {0,0,0,0, "Options to specify when and where on earth" },
  {"location", 'l', "+DDMM+DDDMM|+DDMMSS+DDDMMSS", 0,  "Calculate for latitude (+N/-S) and longitude (+E-W) in Degrees, Minutes and Seconds. Overrides configuration files /etc/zon.conf and ~/.config/zon.conf" },
  {"date",     'd', "YYYY-MM-DDTHH:MM+ZZ:zz", 0,  "Calculate for specified iso-formatted time. Defaults to current system UTC date. Colon in timezone may be omitted (compatibility)" },
  {0,0,0,0, "Options to format the output. Defaults to iso-format" },
  {"at",       '@', 0,      0,  "Format output as date usable by the at command (in UTC), HH:MM YYYY-MM-DD" },
  {"systemd",  'y', 0,      0,  "Format output as required for systemd-run,  YYYY-MM-DD HH:MM UTC" },
  {"format",   'f', "%H %M %m etc",      0,  "Format output yourself with %H:%M %Y-%m%d %Z etcetera, see strftime() documentation" },
  {0,0,0,0, "" },
  {"verbose",  'v', 0,      0,  "Produce a label or if repeated give all base and calculated data, including date and location" },
  {0,0,0,0, "Options to select the kind of twighlight" },
  {"sun",       0,  0,      0,  "Default: Produce start, ending and duration of visibility of top of sun above horizon, i.e. sunrise and sunset" },
  {0,0,0,0, "" },
  {"civil",     1,  0,      0,  "Produce data about civil twighlight, starting when centre of sun is 6 degrees below horizon" },
  {0,0,0,0, "" },
  {"nautical",  2,  0,      0,  "Produce data about nautical twighlight, starting when centre of sun is 12 degrees below horizon" },
  {0,0,0,0, "" },
  {"astronomical",3,0,      0,  "Produce data about astronomical twighlight, starting when centre of sun is 18 degrees below horizon" },
  {0,0,0,0, "" },
  {"angle"     ,5,  "degrees",      0,  "Specify your own rise/set angle of the centre of the sun below horizon" },
  {"rim"       ,4,  0,      0,  "Specify to compensate angle for upper rim of the sun, so to act like sun rise/set" },
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
  regex_t regex;

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
      if (regcomp(&regex, "^((19)|(20))[0-9]{2}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}(\\+|-)[0-9]{2}:?[0-9]{2}$", REG_EXTENDED ) || regexec(&regex, arg, 0, NULL, 0)) argp_usage (state);
      regfree(&regex);
      arguments->date = arg;
      break;
    case 'l':
      if (regcomp(&regex, "^(\\+|-)[0-9]{4}([0-9]{2})?(\\+|-)[0-9]{5}([0-9]{2})?$", REG_EXTENDED ) || regexec(&regex, arg, 0, NULL, 0)) argp_usage (state);
      regfree(&regex);
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
      if ( ! (arguments.current || arguments.rise || arguments.set || arguments.mid ) ) arguments.current=1;

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
	      strncpy(filename,getenv("HOME")?getenv("HOME"):"",80);
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
      time_t tnow,tbase,trise,tset,validtrise,validtset;
      struct tm base,tmrise,tmset;
      char *datestr; datestr=malloc(sizeof(char)*81);

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

      int skipped_days, savedrs;
	   skipped_days = 0 ; 
           validtrise=0;
           validtset=0;
           do {
	     rs = __sunnoonarct__(base.tm_year+1900,base.tm_mon+1,base.tm_mday + skipped_days,lon,lat,  arguments.angle, arguments.rim, &tnoon, &tarc ); 

	     trise=tbase + skipped_days*24*60*60;
	     tmrise= *gmtime(&trise);
	     tmrise.tm_hour = (int)((tnoon-tarc)*60) / 60;
	     tmrise.tm_min  = (int)((tnoon-tarc)*60) % 60;
	     tmrise.tm_sec = 0 ; 
	     tmrise.tm_isdst = 0 ; 
	     trise=timegm(&tmrise); 

	     tset=tbase + skipped_days*24*60*60;
	     tmset= *gmtime(&tset);
	     tmset.tm_hour = (int)((tnoon+tarc)*60) / 60;
	     tmset.tm_min  = (int)((tnoon+tarc)*60) % 60;
	     tmset.tm_sec = 0 ;
	     tmset.tm_isdst = 0 ;
	     tset=timegm(&tmset); 

	     if (skipped_days==0) {
		     savedrs=rs;
		     if (arguments.verbose >= 2) {
	               if ( rs>0 ) printf("++ up entire solar day %+i\n",rs);
	               if ( rs<0 ) printf("-- down entire solar day %+i\n",rs);
	               if ( rs==0 ) printf(" Day with at a sun set and/or rise event. RS %+i\n",rs);
		     }; 
	     };
             if (rs==0) {
	       if (savedrs>0) {
		      if ( (!validtset) && tset>tbase ) validtset=tset;
	              if ( (!validtrise) && validtset && (trise>validtset) ) validtrise=trise;
               } else if (savedrs<0 ){
	              if ( (!validtrise) && trise>tbase ) validtrise=trise;
	              if ( (!validtset) && validtrise && (tset>validtrise) ) validtset=tset;
               } else {
	              if ( (!validtset) && (tset>tbase) ) validtset=tset;
	              if ( (!validtrise) && (trise>tbase) ) validtrise=trise;
	       };
	     };
	 //    printf("RS value : %+i\n",rs);
	 //    printf("skipped_days value : %+i\n",skipped_days);
	 //    printf( "trise  %s", ctime(&trise));
	 //    printf( "tset  %s", ctime(&tset));
	 //    printf( "vtrise  %s", ctime(&validtrise));
	 //    printf( "vtset  %s", ctime(&validtset));
	     skipped_days +=1 ; 
           } while  ( ( (!validtrise)|| (!validtset)) && (skipped_days<365)); 

           tset=validtset;
	   trise=validtrise;
	   
	   if (arguments.rise) {
		   gmtime_r(&trise,&tmrise) ;
		   strftime(datestr,80,arguments.dateformat, &tmrise);
		   printf("%s%s\n",datestr,(arguments.verbose>=1)?" rise":"");
	   }
	   
	   if (arguments.set) {
		   gmtime_r(&tset,&tmset);
		   strftime(datestr,80,arguments.dateformat, &tmset); 
		   printf("%s%s\n",datestr,(arguments.verbose>=1)?" set":"");
	   }

	   if (arguments.current) if (tset<trise)
		   printf("+%s\n",(arguments.verbose >=1)?" up now":"");  
	   else 
		   printf("-%s\n",(arguments.verbose >=1)?" down now":"");  
	   double RAss,decss,rss,azss,altss,d;
	   if (arguments.current) if (arguments.verbose >=2) {
		    d = days_since_2000_Jan_0(base.tm_year+1900,base.tm_mon+1,base.tm_mday) + base.tm_hour/24.0 + base.tm_min/(24*60.0); 
		   sun_RA_dec(d,&RAss,&decss,&rss);
		   EqAz(RAss,decss,base,lon,lat,&azss,&altss);
		   printf("sun azimuth=%f  altitude=%f\n",azss,altss);
	   }

	
	   if (arguments.mid) {
		   tset = (trise+tset)/2;
		   gmtime_r(&tset,&tmset);
		   strftime(datestr,80,arguments.dateformat, &tmset); 
		   printf("%s%s\n",datestr,(arguments.verbose>=1)?" mid":"");
	   }

	//   rs =  __sunrise__( base.tm_year+1900, base.tm_mon+1, base.tm_mday, lon, lat,
        //          arguments.angle, arguments.rim, &tnoon, &tarc );
	//	   printf("trise=%f\n",tnoon);
      return 0;
}

