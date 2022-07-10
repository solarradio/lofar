pro lofar_sunposJ2000, UT ,ra, dec, pa, print=print
; written by Eduard@Glasgow June 2017
; 
; input:
; UT time string e.g. '18:15:00 04-Jun-2015'
; lofar_sunposJ2000, '18:15:00 04-Jun-2015' ,ra, dec, pa, /print
;output
; ra and dec in degrees
; PA angle 
IF (N_elements(PRINT) EQ 0) then print = 0

jd= anytim2jd(UT)
; calulates JD 
;sunpos, jd.int+jd.frac, ra, dec ;& print,ra,dec
; calulates apparent RA and DEC for the sun 
;if (print EQ 1) THEN print, 'Apparent RA =', ra,'  DEC =', dec

;a=jpl_eph(jd.int+jd.frac, 3,center=0)

a=get_sun(UT)

ra=a[8]*15.0D & dec =a[9] & pa=a[12]
;aparrent RA &DEC of the sun
;ra=a[6]*15.0D & dec =a[7] & pa=a[12]
;true RA &DEC of the sun

;ephoch = ( - 2451545.0)/365.25
J = 2000.0d + (JD.int+jd.frac-2451545.0d0)/365.25d
;stop 

precess, ra, dec , J, 2000.0
;sixty(ra)
if (print EQ 1) THEN print, '(J2000) RA =', ra,' degrees  DEC =', dec ,' degrees'
;if (print EQ 1) THEN print, sixty(ra*24/360), FORMAT='The values are: {%f} {%f} {%f}'
if (print EQ 1) THEN print, sixty(ra*24/360), FORMAT = '("J2000 RA =",("", I3, ":"),("", I2, ":"),("", F6.3, " [hh:mm:ss]sec"))'
if (print EQ 1) THEN print, sixty(dec), FORMAT = '("J2000 DEC=",("", I3, " deg "),(" ", I2, " min"),(" ", F6.3, " sec"))'
;stop
end

