
;+
; NAME:
;	lofar_xyz2enu_coordinates
; PURPOSE:
;	Calculates the 3D Cartesian coordinates of a given point at the Earth,
;	as observed from the Sun at a given time.
; CALLING SEQUENCE:
;	ITRF2SolarView, time, ITRF, xyz
; INPUTS:
;	time	reference date and time (in any format acceptable for ANYTIM).
;	ITRF	coordinates of the specified point (points) in meters 
;		according to the International Terrestrial Reference Frame 
;		system (http://itrf.ensg.ign.fr).
;		Either a 3-element array or a [3, N] array where N is the 
;		number of points.
; OUTPUT:
;	xyz	[x, y, z] coordinates of the specified point (points) in 
;		meters at the specified time in the solar-oriented reference 
;		frame, i.e.:
;		reference frame center: the Earth center;
;		z axis: directed to the Sun;
;		y axis: normal to z, directed to the north;
;		x axis: normal to z and y, directed to the east.
;		The output has the same format as the ITRF input parameter.
; ROUTINES CALLED:
;	jd, sunpos, ct2lst
; MODIFICATION HISTORY:
;	14 Sep 2016: Created by Alexey Kuznetsov.
;-

pro lofar_xyz2enu_coordinates, time, ITRF, xyz
 jd=tim2jd(time) ;converts date and time to Julian day number
 sunpos, jd, ra, dec ;calculates the solar right ascension and declination 
                     ;(both in degrees)

 lat=dec ;solar latitude equals its declination

 ct2lst, lst, 0.0, 0.0, jd ;calculates the Local Sidereal Time at Greenwich 
                           ;(in hours)
                           ;LST equals the right ascension of the local zenith
 lon=ra-lst*15 ;solar longitude is the difference between its right ascension 
               ;and the right ascension of the Greenwich zenith

 ;conversion to radians
 lat=lat*!dpi/180
 lon=lon*!dpi/180
 
 ;rotation matrix
 A=[[         -sin(lon),           cos(lon),      0.0], $
    [-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)], $
    [ cos(lat)*cos(lon),  cos(lat)*sin(lon), sin(lat)]]
 
 xyz=transpose(A) # ITRF
 
end