;+
; NAME:
;       lofar_fft_uv2xy
; PURPOSE:
;       calculate dirty map/beam and report information of scales
;       set as a result of uv gridding
; CATEGORY:
;   OVSA APC imaging
; CALLING SEQUENCE:
;      ovsa_fft_uv2xy,f_GHz, m, vis_ij, uv_ij, xyint, dbeam, dmap
; INPUTS:
;       vis_ij : complex visibility as a fcn of time(i) & bsl(j)
;       uv_ij  : complex (u,v) at the corresponding time & bsl
;   m      : dimension of desired map = m x m
; OPTIONAL (KEYWORD) INPUT PARAMETERS:
;       -
; ROUTINES CALLED:
;       -
; OUTPUTS:
;        xyint, dbeam, dmap
;
; COMMENTS:
;    1. Why so named? This code would be similar to fftmap.f
;       ovsa_fft_uv2xy = fft (visibility in uv plane) to xy plane,
;           in clear contrast to
;   ovsa_fft_xy2uv = fft (map in xy plane) to uv plane.
;    2. uv gridding is done here.
;   uvint is automatically set rather than user's choice
;   convol. griid fcn is also fixed to one type: tr. gaussian.
;   
;   For details see here
;   https://casa.nrao.edu/docs/taskref/sdgridold-task.html
;   
;   While dmap is mxm arrary, dbeam is double sized, i.e., 2mx2m
;   to do better cleaning.
; SIDE EFFECTS: -
; RESTRICTIONS: see comments
; MODIFICATION HISTORY:1 by JL
;     15-Dec-2002  DG
;       Changed the sign of the FFT from -1 to the correct +1.
;       All the maps were coming out inverted in x and y.
;     01-Mar-2003  DG
;       Mask uv points that have zero visibilities (otherwise
;       the beam was different than the map when there were
;       missing data).
;     06-Dec-2004  DG
;       Added comment defining the factor 3.6e-5, which is correct!
;   13-March-2007 SDT
;   Changed the normalization scheme for the gridding interpolation scheme.
;   The convolving gaussian is now normalized by its integral and the effect
;   of the convolution is taken into account by dividing the maps by the FT
;     of this normalized gaussian.
;   20-May-2007 SDT
;   Removed the conversion to Tb from this program. The older version did not
;   take into account the division by the beam volume in steradians, but also
;   I wanted to keep all the conversions in one place, as they were previously
;   distributed in several programs.
;     26-Feb-2014 GN
;       Introduced selector keyword
;     14-Mar-2014 GN
;       Renamed to ovsa_ovsa_fft_uv2xy to avoid potential conflict with a copy of the original version
;       located in the ssw\hessi three. All ovsa calling routines were updated accordingly.
;   14-Mar-2014  DG
;      Changed some routine names by appending ovsa_ to avoid potential conflict with the copies of the original
;      routines located in the ssw\hessi tree
;     17 -Feb -2017 Eduard@glasgow
;     Trunc. Gaussian 7 component calculations rewritten to arbitrary number NAA
;     now NAA can be changes without change of the code
;     works good for LOFAR with NAA=23
;-

pro lofar_fft_uv2xy, f_GHz, m, vis_ij, uv_ij, xyint, dbeam, dmap,wsum,skipbls=skipbls, noprint=noprint,selector=selector,vis=vis,uv=uv

  ;Here the user can modify the visibilities used in the clean process.
  ;There are several examples which can be used, but the user must decide
  ;which antennas must be eliminated for the specific data set to be useful.
  ; For gridding purposes, reform ij to 1D array and form
  ; complex conjugate visibilities at mirror-symmetric (u,v)
  ; the matrix indeces are
  ;1  2  4  5  6  7  8
  ;-------------------
  ;-  0  1  2  3  4  5
  ;   -  6  7  8  9 10
  ;      - 11 12 13 14
  ;         - 15 16 17
  ;            - 18 19
  ;               - 20
  ;                  -

  if keyword_set(skipbls) then begin
    goto,top
  endif else begin
    if ~keyword_set(noprint) then print,' '
    if ~keyword_set(noprint) then print,'  FTMAP: UV_gridding and FFT '
    if ~keyword_set(noprint) then print,'  ========================== '
    ;vis_sav = vis_ij
    ;Eliminate Antennas 5,8
    ;vis_ij = vis_ij[*,[0,1,3,4,6,8,9,12,13,15,16,18]]
    ;uv_ij = uv_ij[*,[0,1,3,4,6,8,9,12,13,15,16,18]]
    
    ;Select baselines if selector keyword provided, or select all baselines if selector=-1
    if n_elements(selector) gt 0 then begin
      if selector[0] ne -1 then begin
        vis_ij = vis_ij[*,selector]
        uv_ij = uv_ij[*,selector]
      end
    end
  endelse

  top:
;stop
  nt=n_elements(uv_ij(*,0))
  ; Ask user about removing certain antennas/baselines
  ;junk = ''
  ;read,'Enter antennas (indexes) to remove [e.g. 0,1,4; -1 for none]: ',junk
  ;xblist = strtok(junk,',',/extract)
  xblist = -1
  if (n_elements(xblist) eq 1 and xblist[0] eq -1) then begin
    nb=n_elements(uv_ij(0,*))
    ntb=long(nt)*long(nb)
    uv =reform( uv_ij,ntb) & uv =[ uv,    -uv  ] ; mirror symmetry
    vis=reform(vis_ij,ntb) & vis=[vis,conj(vis)] ; cmplx conjugate
    ;vis=reform(vis_ij,ntb)
    ;vis=reform(vis_ij,ntb) & vis=[vis,0.] ; cmplx conjugate
  endif else begin
    nb=n_elements(uv_ij[0,*])
    blist = indgen(nb)
    blist[xblist] = -1
    blist = blist(where(blist eq -1))
    nb=n_elements(blist)
    ntb=long(nt)*long(nb)
    uv =reform( uv_ij[*,blist],ntb) & uv =[ uv,    -uv  ] ; mirror symmetry
    vis=reform(vis_ij[*,blist],ntb) & vis=[vis,conj(vis)] ; cmplx conjugate
    ;vis=reform(vis_ij[*,blist],ntb)
    ;vis=reform(vis_ij[*,blist],ntb) & vis=[vis,0.]
  endelse

  ;if ~keyword_set(noprint) then print,' Total # of u,v points  :',ntb
  u=float(uv) & v=imaginary(uv)
  tot_uv_siz = 4.0 * max([max(u),max(v)]) ;wanted a dbeam twice as large as needed
  ;to ease shifting without wrapping.
  ; double the mapsize for cleaning in wide enough area
  m2  = m*2
  uvdens  = intarr(m2,m2)
  visarr  = complexarr(m2,m2)
  visbrr  = complexarr(m2,m2)

  ; set xy scale
  ; If XYINT is not zero, use it to determine scales, otherwise
  ; determine scales automatically.
  if (xyint eq 0.) then begin
    xyint = 206265./tot_uv_siz      ; in units of (arcsec)
  endif else begin
    tot_uv_siz = 206265./xyint
  endelse

  uvint = tot_uv_siz / float(m2)

  ;if ~keyword_set(noprint) then print,' Convol. gridding fcn = truncated gaussian'
  NAA=23
  xcon    = fltarr(NAA,NAA)
  for i=0, NAA-1 do xcon(i,*)=float(i-NAA/2)
  ycon    = transpose(xcon)
  fden    = indgen(NAA,NAA)*0+1L

  ;fcon_yn='yes'
  fcon_yn='no'

  case fcon_yn of
    'yes': begin
      if ~keyword_set(noprint) then print,'  Convolution Gridding function  : Trunc. Gaussian'
      ;;is=intarr(n_elements(u)) & js=is
      for ip=0l,n_elements(u)-1 do begin
        xi = u(ip)/uvint+m & i = fix(xi) & di = xi - i
        yj = v(ip)/uvint+m & j = fix(yj) & dj = yj - j
        ;; is[ip]=i & js[ip]=j
        fcon= 1.0 / exp( (xcon-di)^2 + (ycon-dj)^2 )
        fcon=fcon/total(fcon)   ;normalize integral, not peak to unity (SDT)
        visarr(i-NAA/2:i+NAA/2,j-NAA/2:j+NAA/2)= visarr(i-NAA/2:i+NAA/2,j-NAA/2:j+NAA/2) + fcon*vis(ip)
        visbrr(i-NAA/2:i+NAA/2,j-NAA/2:j+NAA/2)= visbrr(i-NAA/2:i+NAA/2,j-NAA/2:j+NAA/2) + fcon*complex(vis[ip] gt 0.,0.)
        uvdens(i-NAA/2:i+NAA/2,j-NAA/2:j+NAA/2)= uvdens(i-NAA/2:i+NAA/2,j-NAA/2:j+NAA/2) + fden*(vis[ip] gt 0.)

      endfor

      ok=where(uvdens gt 0)
      visarr(ok)=visarr(ok)/uvdens(ok)
      visbrr(ok)=visbrr(ok)/uvdens(ok)

      Cfn=fltarr(m2,m2)
      fcon= 1.0 / exp( (xcon)^2 + (ycon)^2 )
      fcon=fcon/total(fcon)
      Cfn[m2/2-NAA/2:m2/2+NAA/2,m2/2-NAA/2:m2/2+NAA/2]=fcon
      FTCfn=shift(fft(shift(Cfn,m2/2,m2/2),1),m2/2,m2/2)

      map  = shift(float(fft(shift(visarr,m,m),1)),m,m)/float(FTCfn)
      beam = shift(float(fft(shift(visbrr,m,m),1)),m,m)/float(FTCfn)
      ;goto,top
    end

    'no': begin

      for ip=0l,n_elements(u)-1 do begin
        i = round(u(ip)/uvint+m)
        j = round(v(ip)/uvint+m)
        visarr(i,j)= visarr(i,j) + vis(ip)
;        stop
        visbrr(i,j)= visbrr(i,j) + complex(1.,0.)
        uvdens(i,j)= uvdens(i,j) + 1L
      endfor

      ok=where(uvdens gt 0)
      visarr(ok)=visarr(ok)/uvdens(ok)
      visbrr(ok)=visbrr(ok)/uvdens(ok)

      map  = shift(float(fft(shift(visarr,m,m),1)),m,m)
      beam = shift(float(fft(shift(visbrr,m,m),1)),m,m)

    end
  endcase
;stop
  ;wsum = beam(m,m)

  ; at this point want to return dbeam and beam are
  ; normalized their max equalt to unity, separating true
  ; value aside
  ; the pkvalu returned is TB_max of the map (in MK units)
  ; except it should be divided by the cleam beam area
  ; in pixel units.
  ; so apply some physical parameters
  ; Note: 3.6e-5 comes from (1.38e-23)*1.e22*1.e6*1.e18/((206265.*3.e8)^2)
  ;                             k      1/SFU  MK  GHz^2   "/rad    c
  ;       This is per polarization.

  wsum=beam(m,m)

  beam = beam/wsum
  map  = map/wsum ; map is very small ~ 10^-4 so mult earlier
  ;;save the conversion to Tb for ovsa_clean_only.pro, for now (SDT)
  ;map  = (map/wsum) /3.6E-5/(xyint*f_GHz)^2
  ;map  = map /(3.6E-5*f_GHz^2)*(uvint/206265.)^2*float(m2)^2; this should be same as the above
  ;if ~keyword_set(noprint) then print,1/(3.6E-5*f_GHz^2)*(uvint/206265.)^2*float(m2)^2,'test tb2'  ;test
  ;map  = map /(3.6E-5*f_GHz^2)*(uvint/206265.)^2*float(m)^2 ;this is 1/tb2flux ,SDT
  ;bmap=map
  map  = map(m/2:3*m/2-1,m/2:3*m/2-1)  ;going down to the smaller map size

  ;flux2tb=1./(3.6E-5*f_GHz^2)*(uvint/206265.)^2*float(m2)^2 ; Need for proper visib retrieval, SDT

  ; To finally get TB in MK we need to divide it by f_GHz^2 and
  ; (pi*a*b) which is done later in the main program (imgr.pro).
  ;UPDATE: The above operation involving 3.6E-5 already divides by f_GHz^2
  ;and it is in ovsa_clean_only.pro that the map is divided by the clean beam, (pi*a*b), SDT
  
  dmap  = map
  dbeam = beam
  xyint = xyint

  if ~keyword_set(noprint) then print,'  Wsum (Sum of weighted gridding): ',wsum
  if ~keyword_set(noprint) then print,'  XYint [arcsec]                 : ',xyint

  return

end
