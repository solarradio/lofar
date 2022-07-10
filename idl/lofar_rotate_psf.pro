pro lofar_rotate_psf
  ; rotates PSF
  ;coordinates of LOFAR core
  lat=52.914610 & lon=6.871362
  
  ; Source *********************
  ;TAU A RA: 83.6333 DEC: 22.0144
   raTauA=83.6333 & decTauA=22.0144
   
   ;time for calculation UT
   obs_time='12-May-17 11:35:00'
   jd0=tim2jd(obs_time)
   eq2hor, ra0, dec0, jd0, alt0, az0,   LAT= lat, LON=lon
   eq2hor, ra1, dec1, jd1, alt1, az1, LAT= lat, LON=lon
   
   xx=(alt0-alt0[0])/cos(!PI*(90-alt[0])/180)*cos(!PI*(90-alt1[0])/180)
   hor2eq, alt, az, jd, ra1, dec1, lat=lat, lon=lon
   
   
  
  
  print, 'ALT=', alt, '     AZ=', az
  
  ;test the reverse
  hor2eq, alt, az, jd, ra1, dec1, lat=lat, lon=lon
  print,RA1-RA,DEC1-DEC1
stop

end
