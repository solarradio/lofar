

pro lofar_radecJ2000_2solar_xy, UT,ra,dec,p0, x,y, tauA=TauA
; 
; example 
; lofar_radecJ2000_2solar_xy, '18:15:00 04-Jun-2015', ra,dec,p0,x, y 
; ra, dec, calculates x, y from the Sun centre
;inputs: ra, dec of a point or array
; x,y coordinates at the Sun for time UT
; inputs are in degrees assuming J2000 coordinates
; created by Eduard@astro for LOFAR project
; http://www.issibern.ch/teams/lofar/
; For details see http://adsabs.harvard.edu/abs/2017NatCo...8.1515K


lofar_sunposJ2000, UT ,raS, decS, p0, /print
; calculates ra, dec, pa for the Sun in J2000 coordinates

; IMPORTANT USING TAU A instead of the SUN
IF N_elements(tauA) EQ 1 THEN BEGIN
raS =83.6331
decS=22.0145
END
;*** should be commented for solar observations

;raS=ra[0]
;decS=dec[0]
;sets the first beam as the Sun 

dec0_rad=decS*!PI/180
ra0_rad = raS*!PI/180


dec_rad=dec*!PI/180
ra_rad = ra*!PI/180

;dec0_rad = (decd + decm/60. + decs/3600.)*!dtor
x_geo = (ra0_rad - ra_rad)*cos(dec0_rad)*206265.80
y_geo = (dec_rad - dec0_rad)*206265.80
p0_rad = p0*!dtor
x = cos(p0_rad)*x_geo + sin(p0_rad)*y_geo
y = -sin(p0_rad)*x_geo + cos(p0_rad)*y_geo

  return
end
