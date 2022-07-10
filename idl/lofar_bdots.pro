;+
; NAME:
;     BDOTS
; PURPOSE:
;     Given the current position coordinates (HA and Dec) and a header
;     containing baseline info, calculate the time delay [nsec].  The
;     routine gets its name from the fact that time delay is tau = B . S
;     where B is the baseline vector and S is the position vector.
; CATEGORY:
;     OVRO SPAN DATA-ANALYSIS
; CALLING SEQUENCE:
;     tau = lofar_bdots(ha,dec,geometry,uv)
; INPUTS:
;     ha        the current hour angle [msec]
;     dec       the current declination [masec]
;     geometry  the geometry structure as returned from NEWSCAN.
; OPTIONAL (KEYWORD) INPUT PARAMETERS:
; ROUTINES CALLED:
; OUTPUTS:
;     tau       an array of size NANT x NANT
;                 containing the delay in nsec for each baseline
;     uv        an NANT x NANT array containing the u,v coordinates in
;                 nsec.  U is contained in the upper non-diagonal elements
;                 and V is in the lower non-diagonal elements.
; COMMENTS:
; SIDE EFFECTS:
; RESTRICTIONS:
; MODIFICATION HISTORY:
;     Written 11-Mar-1997 by Dale Gary
;     10-Jan-1999  DG
;       Changed input argument name from HEADER to GEOMETRY
;     01-Nov-1999  DG
;       Add calculation of u,v coordinates (these will correspond
;       to the start of the cycle.
;     14-Mar-2000  DG
;       Made some things explicitly DOUBLEs, possibly better accuracy.
;-


function lofar_bdots,ha,dec,geometry,uv

  ; Some parameters
  CMAS2RAD = !dtor / 3.6D6
  CMS2RAD = CMAS2RAD * 15.D
  
  ;coordinates of LOFAR core
  lat=52.914610 & lon=6.871362
  CSLAT    = cos(!Pi*lat/180)         ; Cosine of LOFAR latitude 
  SNLAT    = sin(!Pi*lat/180)         ; Sine of LOFAR latitude 
  
  
  ; Convert inputs to radians

  ha1 = ha*CMS2RAD
  dec1 = dec*CMAS2RAD

  ; Determine Cosine of Zenith angle

  cosza = sin(dec1)*SNLAT + cos(dec1)*CSLAT*cos(ha1)
  tau   =   geometry.bx*cos(dec1)*cos(ha1) $
    - geometry.by*cos(dec1)*sin(ha1) $
    + geometry.bz*sin(dec1) + geometry.htdiff/cosza

  ; Calculate the u,v coordinates (this is a 5x5 array, with the upper off-diagonal
  ; elements containing u, and the lower off-diagonal elements containing v)
  uv = geometry.bx*sin(ha1) + geometry.by*cos(ha1)    ; This is just u, but v is added below
  v =   -geometry.bx*sin(dec1)*cos(ha1) $
    + geometry.by*sin(dec1)*sin(ha1) $
    + geometry.bz*cos(dec1)

  ; Get the indexes of the lower off-diagonal elements
  nant = n_elements(uv[0,*])   ; Find out how many antennas from the size of uv
  lower = nant
  for i = 2, nant-1 do begin
    for j = 0, i-1 do begin
      lower = [lower,i*nant+j]
    endfor
  endfor
  uv[lower] = v[lower]                            ; Put v in the lower elements
;stop
  return,tau
end