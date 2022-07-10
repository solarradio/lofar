;+
; NAME:
;     LOFAR_CORE_UV
; PURPOSE:
;     Returns OVSA UV coverage for given HA range, Dec, and antennas
; CATEGORY:
;     OVSA MISC
; CALLING SEQUENCE:
;     uvt = ovsa_uv(harange,hstep,dec,antlist[,halist])
; INPUTS:
;     harange   a 2-element array giving the start and end hour-angle, [h].
;                 Valid OVSA coverage is limited to the range [-4.0,4.0]
;     hstep     a float giving the step size in hour angle [h]
;     dec       the source declination (assumed constant) [deg]
;     antlist   an integer array specifying the OVSA antennas for which to
;                 return coverage, e.g. [1,2,4,5,6,7] for all antennas of
;                 current array.  The 7th antenna's approximate location
;                 can be specified by antenna number 3.
; OPTIONAL (KEYWORD) INPUT PARAMETERS:
; ROUTINES CALLED:
;     uvt       a float array of size [NANT,NANT,NSTEPS], where NANT is
;                 n_elements(antlist), and NSTEPS is derived from the HA
;                 range and step size.  The upper non-diagonal elements
;                 for each step contain the U values, while the lower
;                 non-diagonal elements contain the V values.  The units
;                 can be thought of as nsec, or wavelengths at 1 GHz.
;                 A convenient method for plotting the values for all
;                 returned antennas is:
;                    plot,[-2000,2000]*2,[-2000,2000]*2,/nodata
;                    for i = 0, nant-1 do for j = i+1,nant-1 do $
;                           oplot, uvt[j,i,*], uvt[i,j,*],psym=3
;                    for i = 0, nant-1 do for j = i+1,nant-1 do $
;                           oplot,-uvt[j,i,*],-uvt[i,j,*],psym=3
;     halist    a float array giving the list of hour angles at which the
;                 UV points are determined.  Optional.
; OUTPUTS:
; COMMENTS:
; SIDE EFFECTS:
; RESTRICTIONS:
; MODIFICATION HISTORY:
;     using ovsa_uv by Dale Gary
;     changed 7-Oct-20 by eduard.kontar@glasgow 
;     
;-

function lofar_core_uv,harange,hstep,dec,lon,lat
  
  ;dec=lat
  
  ;conversion to radians
  lat=lat*!pi/180
  lon=lon*!pi/180

  ; Constant parameters
  CMM       = 2.997925D2          ; Speed of light [mm/nsec]
  ;CSLAT    = 0.79619623D         ; Cosine of OVRO latitude (37d13'53.8")
  ;SNLAT    = 0.60503848D         ; Sine of OVRO latitude (37d13'53.8")
  
  ;coordinates of Amsterdam
  ;52.3667° N, 4.9000° E
  ;lat=!pi* 52.3667/180 ;latitude relative to the equator, positive to north, negative to south
  ;lon=!pi*(+4.900)/180 ;longitude relative to the zero meridian, positive to east, negative to west
  
  ;coordinates of LOFAR core 
  ;52.914610, 6.871362

  ;lat=!pi*52.914610/180.0D
  ;lon=!pi*6.871362/180.0D
  
   COSLAT=COS(LAT)
   SINLAT=SIN(LAT)
   COSLON=COS(LON)
   SINLON=SIN(LON)
   ;stop
   ;Nominal XYZ coordinates of LBA stations
   
   CS001= [3826923.546,460915.441,5064643.489]*1d
   CS002= [3826577.066,461022.948,5064892.786]*1d
   CS003= [3826516.748,460930.066,5064946.457]*1d
   CS004= [3826654.197,460939.576,5064842.426]*1d
   CS005= [3826668.750,461069.550,5064819.754]*1d
   CS006= [3826596.730,461145.178,5064866.978]*1d
   CS007= [3826533.361,461098.966,5064918.721]*1d
   CS011= [3826667.069,461285.849,5064801.592]*1d
   CS013= [3826346.265,460792.111,5065087.136]*1d
   CS017= [3826462.054,461501.950,5064935.827]*1d
   CS021= [3826406.543,460538.604,5065064.870]*1d
   CS024= [3827161.234,461409.408,5064421.046]*1d
   CS026= [3826390.916,461869.852,5064955.913]*1d
   CS028= [3825600.445,461260.593,5065604.325]*1d
   CS030= [3826014.266,460387.389,5065372.328]*1d
   CS031= [3826439.996,460273.833,5065063.594]*1d
   CS032= [3826891.573,460387.910,5064715.292]*1d
   CS101= [3825848.343,461689.538,5065378.785]*1d
   CS103= [3826304.279,462823.089,5064934.334]*1d
   CS201= [3826708.929,461913.747,5064713.838]*1d
   CS301= [3827429.462,460990.224,5064257.677]*1d
   CS302= [3827945.916,459792.639,5063990.016]*1d
   CS401= [3826766.106,460100.388,5064836.470]*1d
   CS501= [3825625.779,460642.110,5065640.772]*1d
;     
   
   RPC = [3826577.110,461022.900,5064892.758]*1d
   ; reference phase center
   CS=[[CS001],[CS002],[CS003],[CS004],[CS005],[CS006],[CS007],[CS011],[CS013],[CS017],$
       [CS021],[CS024],[CS026],[CS028],[CS030],[CS031],[CS032],[CS101],[CS103],[CS201],$
       [CS301],[CS302],[CS401],[CS501]]
  
   ; 6-superterp stations
   
   ;  CS=[[CS001],[CS002],[CS003],[CS004],[CS005],[CS006],[CS007],[CS011],[CS013],[CS017],$
   ;    [CS021],[CS024],[CS026],[CS028],[CS030],[CS031],[CS032],[CS101],[CS103],[CS201],$
   ;    [CS301],[CS302],[CS401],[CS501]]
       
   CSNAME=['CS001','CS002','CS003','CS004','CS005','CS006','CS007','CS011','CS013','CS017',$
       'CS021','CS024','CS026','CS028','CS030','CS031','CS032','CS101','CS103','CS201',$
       'CS301','CS302','CS401','CS501']
   ;CSNAME=INTARR(24)
   
   E=dblarr(N_elements(CS(0,*)))
   N=dblarr(N_elements(CS(0,*)))
   U=dblarr(N_elements(CS(0,*)))
   
   ii=[         -sin(lon),           cos(lon),      0.0]
   jj=[-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)]
   kk=[ cos(lat)*cos(lon),  cos(lat)*sin(lon), sin(lat)]
   
   POSITION=CS-REBIN(RPC,3,N_elements(CS[0,*]))
   ; in mm
 
     E=position[0,*]*ii[0]+position[1,*]*ii[1]+position[2,*]*ii[2]
     N=position[0,*]*jj[0]+position[1,*]*jj[1]+position[2,*]*jj[2]
     U=position[0,*]*kk[0]+position[1,*]*kk[1]+position[2,*]*kk[2]
    
    ;stop
    for i=0, n_elements(E)-1 DO BEGIN
      ENU=XYZ2ENU(CS[*,i],RPC)
      E[i]=enu[0] 
      N[i]=enu[1] 
      U[i]=enu[2]
    endfor
  
;  stop
   Antxyz=[E,N,U]*1000.D0
   ; converts to mm
  ; Antxyz=Position*1000.d0
   
   window,xsize=800,ysize=800
   
   plot,E,N, xrange=[-2000,2000],yrange=[-2000,2000],psym=7, xtitle='W-E, [meters]',ytitle='S-N, [meters]'
   plots,[0,0],[-2000,2000],line=1
   plots,[-2000,2000],[0,0],line=1
   x=findgen(2000)-1000+0.5
   oplot,x,sqrt(1e3^2-x^2),line=1
   oplot,x,-sqrt(1e3^2-x^2),line=1
   XYOUTS, E, N,string(CSNAME), CHARSIZE = 1, /data
   ;XYOUTS,200.,200.,'T',CHARSIZE = 1,/data

;  stop
   
  ; Nominal Antenna Positions [mm], in E, N, U coordinates
  ;
  ;     E        N        U
  ; -------- -------- --------
  ;    60960        0        0
  ;   487680        0        0
  ;  -184090     7320   -12680
  ;    13870   128120   -12530
  ;    13870   371960   -12340
  ;   149962   124920   -12500
  ;  1066800   -91440   -12500     ; Approximate location of Ant 3
  ;       -1       -1       -1

;  antxyz = [[  60960,       0,        0],$
;    [ 487680,       0,        0],$
;    [-184090,    7320,   -12680],$
;    [  13870,  128120,   -12530],$
;    [  13870,  371960,   -12340],$
;    [ 149962,  124920,   -12500],$
;    [1066800,  -91440,   -12500],$   ; Approximate location of Ant 3
;    [     -1,      -1,       -1]]

   nant=N_elements(CS[0,*])
   ; all antenas are used
   iant=indgen(N_elements(CS[0,*]))
   ; antenna numbers
   
   
  ; Create GEOMETRY structure needed for BDOTS argument
  geometry = {bx:fltarr(nant,nant),by:fltarr(nant,nant),$
    bz:fltarr(nant,nant),htdiff:fltarr(nant,nant)}

  ; Get current baseline corrections, as [X,Y,Z] in nsec for each of eight
  ; antennas.  Unused antennas have zeroes.  BLCOR is of size [3,8].
  ;blcor = get_blcor()
  
  bcor=0.0
  ; Create matrix of baselines, with convention that baseline
  ; Bx = ANTjx - ANTix. i.e. the baseline vector for baseline ij
  ; goes from ANTi to ANTj.

  for i = 0, nant-1 do begin
    for j = i+1, nant-1 do begin

      ; First get East, North, Up coordinates for baseline ij [nsec]

      be = (antxyz[0,iant[j]] - antxyz[0,iant[i]])/CMM
      bn = (antxyz[1,iant[j]] - antxyz[1,iant[i]])/CMM
      bu = (antxyz[2,iant[j]] - antxyz[2,iant[i]])/CMM

      ; The height difference for each baseline can be used to
      ; correct for atmospheric refraction (refractive index n,
      ; where n - 1 = 0.000285).

      sgn = 1.0
      ;sgn1=0.0
      ;removing the correction
      geometry.htdiff[i,j] = sgn*(bu*0.000285D)
      geometry.htdiff[j,i] = geometry.htdiff[i,j]

      ; Convert baseline coordinates from east, north, up to x,y,z
      ; system of coordinates [nsec]

      bx =   bn*SINLAT - bu*COSLAT
      by = - be
      bz = - bn*COSLAT - bu*SINLAT

      ; Add baseline corrections obtained earlier from BLCIN function.

;      geometry.bx[i,j] = sgn*(bx+(blcor[0,iant[i]]-blcor[0,iant[j]]))
;      geometry.by[i,j] = sgn*(by+(blcor[1,iant[i]]-blcor[1,iant[j]]))
;      geometry.bz[i,j] = sgn*(bz+(blcor[2,iant[i]]-blcor[2,iant[j]]))
     
     geometry.bx[i,j] = sgn*bx
     geometry.by[i,j] = sgn*by
     geometry.bz[i,j] = sgn*bz

      geometry.bx[j,i] = geometry.bx[i,j]
      geometry.by[j,i] = geometry.by[i,j]
      geometry.bz[j,i] = geometry.bz[i,j]
    endfor
  endfor


  ; Determine the number of HA steps to take
  nsteps = fix((harange[1] - harange[0])/hstep)
  if ((harange[0] + nsteps*hstep) le harange[1]) then nsteps = nsteps+1
   ;stop
  ; Declare storage for output
  uvt = fltarr(nant,nant,nsteps)

  ; Convert DEC to msec
  decms = dec*3.6D6
  halist = fltarr(nsteps)

  ; Now loop over HA range, with HASTEP
  for i = 0, nsteps-1 do begin
    ; Determine new HA in hours, and save in HALIST
    ha = float(harange[0]) + i*hstep
    halist[i] = ha
    ; Convert HA to msec
    ha = ha*3.6D6
    ; Call BDOTS repeatedly with HA [msec] and DEC [msec]
    tau = lofar_bdots(ha,decms,geometry,uv)
    uvt[*,*,i] = uv
  endfor
;stop
  return,uvt
end