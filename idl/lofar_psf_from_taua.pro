  pro lofar_psf_from_tauA,obs_time=obs_time,mhz=mhz
  
  If N_ELEMENTS(mhz) LT 1 THEN mhz=31.0
  IF ((MHz GT 50) OR (MHz LE 29)) THEN message, 'Frequency should be 29<MHz<50'
  If N_ELEMENTS(obs_time) LT 1 THEN obs_time='15-Apr-15 11:55:29.821'
  ; rotates PSF
  ;coordinates of LOFAR core
  lat=52.914610 & lon=6.871362
  
  ; Source *********************
  ;TAU A RA: 83.6333 DEC: 22.0144
   raTauA=83.6333 & decTauA=22.0144
   
   ;time for calculation UT
   cd,'D:/lofar'
   restore,filename='psf_tauA_217_beams.sav',/v
   print,'data taken: ',atime(cube_psf.time)
   jd0=tim2jd(cube_psf.time)
   RA0=cube_psf.ra
   DEC0=cube_psf.dec
   
   jd1=tim2jd(obs_time)
   RA1=raTauA & DEC1=decTauA
   
   eq2hor, ra0, dec0, jd0, alt0, az0, LAT= lat, LON=lon
   eq2hor, ra1, dec1, jd1, alt1, az1, LAT= lat, LON=lon
   
   ;xx0=sin((90-alt0)*!PI/180)*cos(az0*!PI/180)
   ;yy0=sin((90-alt0)*!PI/180)*sin(az0*!PI/180)
   ;zz0=cos((90-alt0)*!PI/180)
   
   xx0=(alt0-alt0[0])/cos(!PI*(90-alt0[0])/180)
   yy0=(az0  -az0[0])
   xx=xx0*cos(!PI*(az1[0]-az0[0])/180.)-yy0*sin(!PI*(az1[0]-az0[0])/180.)
   yy=xx0*sin(!PI*(az1[0]-az0[0])/180.)+yy0*cos(!PI*(az1[0]-az0[0])/180.)
   
   alt_P=xx*cos(!PI*(90-alt1[0])/180)+alt1[0]
   az_P =yy+az1[0]
   
   hor2eq, alt_P, az_p, jd1, ra_P, dec_P, lat=lat, lon=lon
   window,xsize=400,ysize=400
   plot,ra_p-ra_p[0],dec_p-dec_p[0],xrange=[-1,1],yrange=[-1,1]
   
   minFreq=min(abs(cube_psf.FREQS-MHZ), iMHZ)
   psf_beams=cube_psf.cube[*,iMHz]
   xout=findgen(200)/100-1.
   yout=xout
   dx=3600*(xout[1]-xout[0])
   dy=dx
   psf_beams=psf_beams-average(psf_beams)
   psf_beams=psf_beams-min(psf_beams)
   psf_beams=psf_beams/max(psf_beams)
   mapxy=GRIDDATA(ra_p-ra_p[0],dec_p-dec_p[0],psf_beams,/RADIAL_BASIS_FUNCTION,xout=xout,yout=yout,VARIOGRAM=[1,8,0,2],/grid)
   mapBeams=make_map(mapxy/max(mapxy), xc=0.,yc=0.,dx=dx,dy=dy,id='LOFAR', SOHO=0B,L0=0,B0=0,RSUN=960)
   plot_map,mapBeams,levels=[0.5,0.7,0.9]*max(mapBeams.data),thick=1,/limb,/bar

  stop
end
