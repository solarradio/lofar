function lofar_sun_zenith_angle, obs_time=obs_time
  ; returns Sun Zenith angle in degrees for LOFAR observations for UT time
  ; INPUT:
  ; obs_time - UT time of the Sun observations 
  
  ; uses SSW routine get_zenang
  
  ; Example:
  ;IDL> lofar_sun_zenith_angle(obs_time='16-Apr-2015, 11:57:00')
  ;
  ; Version 1
  ; Created by Eduard@Glasgow December 20, 2018 
  
   ;time for calculation UT
   IF (N_elements(obs_time) NE 1) then obs_time='16-April-2015 11:57:00'
   
   ;coordinates of LOFAR core
   lat=52.914610 & lon=6.871362
  
   degrees=get_zenang(obs_time, lat, lon, /degrees)
   
   print, ' Zenith_angle=',degrees, ' at UT time =', atime(obs_time)
   return, degrees
   
end
