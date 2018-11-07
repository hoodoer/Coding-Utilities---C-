#ifndef GEOCOORDS_H
#define GEOCOORDS_H

#include <math.h>



/***************************************************************
 Simple include file for handling sphere and spheroid (elliptical)
 conversions between X/Y (in meters) and latitude longitude

 Longitude converts to X (east/west)
 Latitude converts to Y (north/south)

 Don't forget to add -lm to compile options

 Position (0, 0) in XY is lat 0, lon 0.
 Latitude accuracy falls off above 45 degrees and below -45 degrees
*************************************************************/



#define D_R (M_PI / 180.0)
#define R_D (180.0 / M_PI)
#define R_MAJOR 6378137.0
#define R_MINOR 6356752.3142
#define RATIO (R_MINOR/R_MAJOR)
#define ECCENT (sqrt(1.0 - (RATIO * RATIO)))
#define COM (0.5 * ECCENT)



// Max Y's are overly large since
// the equations fall apart above
// 45 degrees lat, and below -45 degrees lat
#define MAX_GEO_X  25000000.0
#define MIN_GEO_X -25000000.0
#define MAX_GEO_Y  50000000.0
#define MIN_GEO_Y -50000000.0







/*******************************************************
// Support functions....
*******************************************************/
inline double deg_rad(double ang)
{
  return ang * D_R;
}


inline double rad_deg(double ang)
{
  return ang * R_D;
}





/*******************************************************
// Earth as sphere conversions
*******************************************************/
inline double sphereXtoLong(double x)
{
  return x / 111120.0;
}


inline double sphereLongToX(double lon)
{
  return lon * 111120.0;

}


inline double sphereYtoLat(double y)
{
  return (R_D * (2 * atan(exp((y / 111120.0) * D_R)) - M_PI / 2));
}


inline double sphereLatToY(double lat)
{
  return 111120.0 * (R_D * log(tan(M_PI / 4 + lat * D_R / 2)));
}


inline void sphereXyToLatLong(double x, double y, double &lat, double &lon)
{
  lat = sphereYtoLat(y);
  lon = sphereXtoLong(x);
}


inline void sphereLatLongToXy(double lat, double lon, double &x, double &y)
{
  x = sphereLongToX(lon);
  y = sphereLatToY(lat);
}






/*******************************************************
// Earth as elliptical spheroid conversions
*******************************************************/
inline double spheroidLongToX(double lon)
{
  return R_MAJOR * deg_rad(lon);
}




inline double spheroidXtoLong(double x)
{
  return rad_deg(x) / R_MAJOR;
}



inline double spheroidLatToY(double lat)
{
  lat = fmin (89.5, fmax (lat, -89.5));
  double phi = deg_rad(lat);
  double sinphi = sin(phi);
  double con = ECCENT * sinphi;
  con = pow((1.0 - con) / (1.0 + con), COM);
  double ts = tan(0.5 * (M_PI * 0.5 - phi)) / con;

  return 0 - R_MAJOR * log(ts); 
}



inline double spheroidYtoLat(double y)
{
  double ts = exp ( -y / R_MAJOR);
  double phi = M_PI_2 - 2 * atan(ts);
  double dphi = 1.0;
  int i;

  for (i = 0; fabs(dphi) > 0.000000001 && i < 15; i++) 
    {
      double con = ECCENT * sin (phi);
      dphi = M_PI_2 - 2 * atan (ts * pow((1.0 - con) / (1.0 + con), COM)) - phi;
      phi += dphi;
    }

  return rad_deg (phi); 
}



inline void spheroidXyToLatLong(double x, double y, double &lat, double &lon)
{
  lat = spheroidYtoLat(y);
  lon = spheroidXtoLong(x);
}


inline void spheroidLatLongToXy(double lat, double lon, double &x, double &y)
{
  x = spheroidLongToX(lon);
  y = spheroidLatToY(lat);
}



#endif
