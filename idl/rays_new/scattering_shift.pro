
;+
; NAME       :  scattering_shift()
; 
; PURPOSE    :  Estimates the radial shift, caused by radio-wave scattering, of the observed radio 
;                 source away from its true location (and away from the Sun), using the simple, 
;                 analytical method derived by Chrysaphi et al. (2018).
;                 
; WRITTEN BY :  Nicolina Chrysaphi @ Glasgow, UK
;
;
; CALLING SEQUENCE :
;       result = scattering_shift(f=f, e2h=e2h [, /harmonic])
;       
; INPUTS:
;       f         --- The observed frequency given in MHz.  A 1D array of frequencies is also accepted.
;       e2h       --- The ratio of epsilon^2/h given in km^(-1).
;                     NOTE: The minimum and maximum values used by Chrysaphi et al. (2018) for 
;                           epsilon^2/h were 7e-5 and 4.5e-5 km^(-1), respectively.
;       
; OPTIONAL INPUT (KEYWORD) PARAMETERS :
;       harmonic  --- If set, harmonic emission is assumed (i.e. f_observed = 2*f_pe).
;                     The default is harmonic=0, meaning that fundamental emission is assumed (i.e. f_observed = f_pe).
;                     
; OUTPUT:
;       Returns the radial shift of the observed radio source away from the true source, in solar radii.
;       Note that scattering in the solar corona shifts sources away from the Sun (and the shift is frequency-dependent).
;       The output has the same data type and dimensions as input f.
; 
; 
; OTHER NOTES :
;       * The analytical method applied assumes that epsilon^2/h is constant with heliocentric distance, and 
;           takes the 1xNewkirk (1961) coronal density model.
;       * This method is valid for all radio emissions, irrespective of their exciter.
;       * For a comprehensive explanation of the method and the assumptions made, 
;           see section 3.4 in Chrysaphi et al. (2018).
;       * Access to the publication by Chrysaphi et al. (2018):
;           Link to the publisher's website:  https://iopscience.iop.org/article/10.3847/1538-4357/aae9e5
;           Link to NASA ADS               :  https://ui.adsabs.harvard.edu/abs/2018ApJ...868...79C
;           Link to arXiv                  :  https://arxiv.org/abs/1810.08026
;       * The scattering_shift() function calls the density_newkirk() function which is included in this routine.
;       * This routine was written in IDL (version 8.3.0).
;  
; MODIFICATION HISTORY :
;       31-Oct-2019 - Nicolina Chrysaphi - Written as an IDL function.
;       14-Nov-2019 - Nicolina Chrysaphi - 1. Adjusted to accept an array of frequencies.
;                                          2. Added COMPILE_OPT IDL2 to all functions.
;       04-Jun-2020 - Nicolina Chrysaphi - Simplified by removing the R_newkirk() function.
;-

FUNCTION density_newkirk, r, N
  ;1xNewkirk (1961) coronal density model:
  COMPILE_OPT IDL2
  n_0 = 4.2e4
  return, N*n_0*10.^(4.32/r)
  ;density in cm^-3
END

FUNCTION scattering_shift, f=f, e2h=e2h, harmonic=harmonic
  
  COMPILE_OPT IDL2
  
  NR = 29000L
  r = findgen(Nr)/Nr*50+1.0
  dr = r[1]-r[0]
  dens = density_newkirk(r, 1.0) ;cm^-3
  fpe = 8.9e-3*sqrt(dens) ;MHz
  
  f = f ;MHz
  eps2oh = e2h  ;km^-1
  Rsun = 6.96d5 ;km
  
  scatt_shift = make_array(n_elements(f),/DOUBLE)
  
  FOR i=0,n_elements(f)-1 DO BEGIN
    fx = fpe[where(fpe LT f[i])]
    rx = r[where(fpe LT f[i])]

    IF keyword_set(harmonic) EQ 1 THEN BEGIN
      tau_dr = (sqrt(!PI)/2.)*eps2oh*(fx^4/(( (2.*f[i])^2 - fx^2 )^2))
      ;For Harmonic emissions: f = 2*fx
    ENDIF ELSE BEGIN
      tau_dr = (sqrt(!PI)/2.)*eps2oh*(fx^4/((f[i]^2-fx^2)^2))
      ;For Fundamental emissions: f = fx
    ENDELSE
    tau = fltarr(n_elements(fx))

    F_GT_FPE = WHERE(fx LT f[i])
    FOR j=min(F_GT_FPE),N_elements(fx)-1 DO BEGIN
      tau[j] = total(tau_dr[j:N_elements(rx)-1])*dr*Rsun
      tau[N_elements(rx)-1] = tau[N_elements(rx)-2]
    ENDFOR

    tau1 = value_locate(tau, 1.0)
    r1 = rx[tau1]

    ;Interpolate between the two closest tau=1 values to find exactly where tau=1:
    INTEP, tau, rx, 1.0, r1_exact

    scatt_shift[i] = ABS(r1_exact - min(rx))
    ;print, 'The radial shift away from the true source (of frequency ', f[i],' MHz) due to scattering effects was estimated to be = ', scatt_shift[i], ' R_solar .............'
  ENDFOR
  
  return, scatt_shift

END

