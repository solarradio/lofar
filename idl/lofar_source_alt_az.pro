pro lofar_source_alt_az
  ; horizon (alt,az) coordinates from equatorial (ra,dec) coords.


   ;time for calculation UT
   obs_time='12-May-17 11:35:00'
   jd=tim2jd(obs_time)
  
  ; Source ********************* 
  ;TAU A RA: 83.6333 DEC: 22.0144
  ;ra=83.6333 & dec=22.0144
  
  ; Sun position RA/DEC J2000
  sunposJ2000, obs_time ,raS, decS, p0, /print
  ra=RaS*[1,1.1] & dec=DECS*[1,1.1]
   
  ;coordinates of LOFAR core
  lat=52.914610 & lon=6.871362
  

  eq2hor, ra, dec, jd, alt, az, LAT= lat, LON=lon
  
  print, 'ALT=', alt, '     AZ=', az
  
  ;test the reverse
  hor2eq, alt, az, jd, ra1, dec1, lat=lat, lon=lon
  print,RA1-RA,DEC1-DEC1
stop

end
