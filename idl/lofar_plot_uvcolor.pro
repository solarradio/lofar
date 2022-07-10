pro lofar_plot_uvcolor,uvt,_extra=_extra

;plots LOFAR UV in colors
;main part taken from OVSA
 
   nant = n_elements(uvt[*,0,0])
   ntimes = n_elements(uvt[0,0,*])
   uls = (vls = fltarr(2*(nant-2),ntimes))
   uss = (vss = fltarr((nant-3)*(nant-2)/2,ntimes))
   k = -1
   k2 = -1
   for i = 0, nant-2 do begin
      for j = i+1,nant-1 do begin
         if (i eq 0 and j eq 1) then begin
            vll = reform(uvt[i,j,*])
            ull = reform(uvt[j,i,*])
         endif else if (i le 1 and j gt 1) then begin
            k = k + 1
            vls[k,*] = reform(uvt[i,j,*])
            uls[k,*] = reform(uvt[j,i,*])
         endif else begin
            k2 = k2 + 1
            vss[k2,*] = reform(uvt[i,j,*])
            uss[k2,*] = reform(uvt[j,i,*])
         endelse
      endfor
   endfor

   uvmax = max(uvt)>(-min(uvt))
   ps = 6
   window,7,xsize=600,ysize=600
   plot,xran=[-uvmax,uvmax],yran=[-uvmax,uvmax],ull,vll,/nodata,psym=ps,xsty=1,ysty=1,_extra=_extra
   oplot,uss,vss,color=250,psym=ps,symsize=0.25
   oplot,-uss,-vss,color=250,psym=ps,symsize=0.25
   oplot,uls,vls,color=210,psym=ps,symsize=0.25
   oplot,-uls,-vls,color=210,psym=ps,symsize=0.25
   oplot,ull,vll,color=50,psym=ps,symsize=0.25
   oplot,-ull,-vll,color=50,psym=ps,symsize=0.25

return
end