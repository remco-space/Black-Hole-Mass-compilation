pro readBHinformation,input,name,type,select,mbh,dmbh,upperlimit,sig,dsig, $
rekpc ,lk ,c28, drekpc ,dlk ,dc28,covar 

BHs=mrdfits(input,1) 
for i=0,N_ELEMENTS(bhs.name)-1 do bhs[i].   name=strn(bhs[i].   name)
for i=0,N_ELEMENTS(bhs.name)-1 do bhs[i].hetname=strn(bhs[i].hetname)
for i=0,N_ELEMENTS(bhs.type)-1 do bhs[i].type   =strn(bhs[i].type)

;help,/struct,bhs
;** Structure <380b208>, 29 tags, length=168, data length=162, refs=1:
;   NAME            STRING    'MW'
;   HETNAME         STRING    'MW'            HETMGS name (vdB+15)
;   SELECTED        INT             -3
;   RA              FLOAT          -1234.00
;   DEC             FLOAT          -1234.00
;   DDIST           FLOAT        0.00800000
;   DISTERR         FLOAT        0.00100000
;   MBH             FLOAT           6.63300
;   DMBH            FLOAT         0.0480000
;   UPPERLIMIT      FLOAT           0.00000   Flag to indicate an BH mass upperlimit
;   SIG             FLOAT           2.00450
;   DSIG            FLOAT         0.0854000
;   LK              FLOAT          -1234.00
;   DLK             FLOAT              -NaN
;   R50             FLOAT          -1234.00
;   DR50            FLOAT              -NaN
;   C28             FLOAT          -1234.00
;   DC28            FLOAT              -NaN
;   LAGN            FLOAT          -1234.00
;   DLAGN           FLOAT              -NaN 
;   LKOR            FLOAT              -NaN  Convience unit of (LK/R50) to remove covariance
;   DLKOR           FLOAT              -NaN
;  Covariance terms:                  
;   CR50LK          FLOAT          -1234.00
;   CR50C28         FLOAT          -1234.00
;   CLKC28          FLOAT          -1234.00
;   CLKORC28        FLOAT          -1234.00
; Housekeeping
;   TYPE            STRING    'star'
;   REF             STRING    'ghez08,gillessen09              '
; HETMGS sigma_c measurements from vdB+15
;   HETMGSSIG       FLOAT          -1234.00
;   HETMGSDSIG      FLOAT          -1234.00

;   SELECTED        INT           1: in fit 0: omitted -3: missing data (table 3) -4: Globular clusters
s=where(bhs.selected ge 0) 
bhs=bhs[s]   	

name       = bhs.name
type       = bhs.type
select     = bhs.selected
mbh        = bhs.mbh
dmbh       = bhs.dmbh
upperlimit = bhs.upperlimit
sig        = bhs.sig
dsig       = bhs.dsig
rekpc      = bhs.r50
lk         = bhs.lk
c28        = bhs.c28
drekpc     = bhs.dr50
dlk        = bhs.dlk
dc28       = bhs.dc28
covar      = bhs.cr50lk  
	
end

pro set_plotstyle
	set_plot,'ps'
	device,/helvetica,/cmyk
	!p.charsize=0.8
	device,xsize=8.7,ysize=8.7
	device,/encap
	!p.thick=2
	!y.margin=[3.5,0.8]
	!x.margin=[6.0,2.0]
end


pro msigma_fit,input

cleanplot,/silent
set_plotstyle 
!y.margin=[3.5,2.0]
!p.multi=[0,1,1]
!p.font=0
!p.charsize=4.0
!p.thick=1.0
set_plot,'ps'
device,/helvetica
!p.charsize=0.8
device,xsize=8.7,ysize=4.7
!p.thick=3
!y.margin=[3.5,0.8]
!x.margin=[5.5,2.0]
	

input='BHcompilation.fits'
readBHinformation,input,name,type,select,mbh,dmbh,upperlimit,sig,dsig, $
rekpc ,mk ,c28, drekpc ,dmk ,dc28,covar

n=N_ELEMENTS(name)



; set vector for the upper-limits/non-dectections
delta=select*0+1
t=where(upperlimit eq 1)
delta[t]=0


s=where(select gt 0,ngal)

; set up the color vector
color=replicate(cgcolor('grey'),n)
t=where(type eq 'maser',c)
if c gt 0 then color[t]=cgcolor('green')
t=where(type eq 'gas' or type eq 'CO',c)
if c gt 0 then color[t]=cgcolor('cyan')
t=where(type eq 'star' or type eq 'stars',c)
if c gt 0 then color[t]=cgcolor('red')
t=where(type eq 'reverb',c)
if c gt 0 then color[t]=cgcolor('black')
t=where(type eq 'stis',c)
if c gt 0 then color[t]=cgcolor('cyan')
t=where(select eq 0,c)
if c gt 0 then color[t]=cgcolor('darkgrey')


; axis ranges
bhrange=[2.0,10.99]
sigmarange=[1.0,2.8]
bhfprange=[8.3,11.3]

;Figure 1: M-sigma
bh1Dlinmix,sig ,dsig,mbh,dmbh,delta,s,name,color,file='fig/bhsigma.eps', $
	xtitle='velocity dispersion (log km/s )',ytitle='BH mass (log M'+sunsymbol()+' )' ,xrange=sigmarange,/xstyle,yrange=bhrange,/ystyle,yunits='log M_{BH} ',xunits='log \sigma '

;Figure 6: The fundamental plane      
bh2Dlinmix,[[mk],[rekpc]] ,[[dmk],[drekpc]] ,sig,dsig,delta*0+1,s,name,color,$
xunits=['log L_k/L'+sunsymbol()+' ','log R_e/kpc '],ytitle='velocity dispersion (log km/s )' ,file='fig/fp.eps', $
yrange=sigmarange,/ystyle,xrange=bhfprange,/xstyle,yunits='log \sigma_e ',post=postfp

;Figure 7: The BH-size-luminosity relation 
bh2Dlinmix,[[mk],[rekpc]] ,[[dmk],[drekpc]] ,mbh,dmbh,delta,s,name,color,covar=covar,$
      xunits=['log L_k/L'+sunsymbol()+' ','log R_e/kpc '],ytitle='BH mass (log M'+sunsymbol()+' )', $
      /ystyle,file='fig/bhlumsize.eps' ,yrange=bhrange ,xrange=bhfprange,/xstyle,yunits='log M_{BH} '


; assume that 3/4 of the tilt is due to ML change. S5.2.1
mlgamma=0.75*mean(2-1/postfp.beta[0]) ; FP tilt
mlalpha=1.0 & mlbeta=166.0 ; ML_k=1.0 @ 166 km/s  
mlcst= alog10(mlalpha) - mlgamma *alog10(mlbeta)  ; alog(M/L)=mlcst + mlgamma*alog(sigma)
ml=(((10d^sig)/mlbeta)^mlgamma)*mlalpha 
totmass=alog10(ml * 10^mk) ; total stellar mass of the galaxies

; BH-size-mass relation equation (5)
bh2Dlinmix,[[totmass],[rekpc]] ,[[dmk],[drekpc]] ,mbh,dmbh,delta,s,name,color,covar=covar,$
      xunits=['log M_{star} ','log R_e '],ytitle='BH mass (log M'+sunsymbol()+' )', $
      /ystyle,file='fig/bhmasssize.eps' ,yrange=bhrange ,xrange=bhfprange,/xstyle,yunits='log M_{BH} '

; fig 2: BH -- Total mass 
bh1Dlinmix,[totmass] ,[dmk],mbh,dmbh,delta,s,name,color,file='fig/bhtotmass.eps', $
xtitle='stellar mass (log M'+sunsymbol()+' )',ytitle='BH mass (log M'+sunsymbol()+' )',yunits='log M_{BH} ' ,xrange=[8.0,12.5],/xstyle,yrange=bhrange,/ystyle,overplotbhlum=1

; BH - total luminosity
bh1Dlinmix,[mk] ,[dmk],mbh,dmbh,delta,s,name,color,file='fig/bhlum.eps', $
xtitle='total luminosity (log(L_k))',ytitle='black hole mass (log(M_sun))' ,xrange=[8.0,12.5],/xstyle,yrange=bhrange,/ystyle,overplotbhlum=1



end




;=========================================
pro bh2Dlinmix,xm,xsig_in,ym_in,ysig,delta,s,htname,post=post,color,covar=covar,file=file,xunits=xunits,_EXTRA=_EXTRA,yunits=yunits

	if not keyword_set(yunits) then yunits=''

nv=2
xsig=xsig_in[s,*]
yvar=ysig[s]^2
ym=ym_in[s]
n=N_ELEMENTS(s)
xvar = dblarr(n,nv,nv) ; nx,np,np with p=2 (x,y)
for j=0,n-1 do begin ; assumes no covariance sig_xy=sig_yx=0
  for i=0,nv-1 do begin
	xvar[j,i,i] = xsig[j,i]^2 ; sig_xx = sigx
endfor
if keyword_set(covar ) then begin
   xvar[j,0,1] = covar[s[j]]
   xvar[j,1,0] = covar[s[j]]
endif
endfor

mlinmix_err, xm[s,*], ym, post, XVAR=xvar, YVAR=yvar, delta=delta[s],/silent ;,NGAUSS=10

alpha=mean(post.alpha)
beta=[mean(post.beta[0]),mean(post.beta[1])]
epsilon=mean(sqrt(post.sigsqr))
varepsilon=stddev(sqrt(post.sigsqr))

bhpred= (BETA/beta[0]) ## xm ;+ EPSILON

t=where(xm gt -1000 and ysig gt -1000 and ysig lt 4 and delta eq 1,nc);  and color gt -1,nc)

if keyword_set(file) then begin
    set_plotstyle
	device,file=file;,ysize=10.0
endif

plot,bhpred[t],ym_in[t],/yno,psym=1,/nodata, $ 
xtitle=textoidl(xunits[0]+strn(beta[1]/beta[0],format='(F6.2)')+xunits[1]) ,_EXTRA=_EXTRA

; plot the regression and the scatter.
oplot,!x.crange,alpha + beta[0]*!x.crange                        ,color=cgcolor('black')
oplot,!x.crange,alpha + beta[0]*!x.crange+epsilon,linestyle=2    ,color=cgcolor('grey')
oplot,!x.crange,alpha + beta[0]*!x.crange-epsilon,linestyle=2    ,color=cgcolor('grey')
oplot,!x.crange,alpha + beta[0]*!x.crange+epsilon*3,linestyle=1  ,color=cgcolor('grey')
oplot,!x.crange,alpha + beta[0]*!x.crange-epsilon*3,linestyle=1  ,color=cgcolor('grey')

if keyword_set(covar ) then  $
   xerr= (sqrt((beta[0]*xsig_in[*,0])^2+(beta[1]*xsig_in[*,1])^2 +2*beta[0]*beta[1]*covar) )/beta[0] $
else $
   xerr=sqrt((xsig_in[*,0])^2+(beta[1]/beta[0]*xsig_in[*,1])^2) 
;https://en.wikipedia.org/wiki/Propagation_of_uncertainty


pointsize=0.4

; find out which error bars to plot
relxerr=xerr/abs(!x.crange[1]-!x.crange[0])
relyerr=ysig/abs(!y.crange[1]-!y.crange[0])
te=where(xm gt -1000 and ysig gt -1000 and  color ne cgcolor('darkgrey'),nc)
srx=te[((sort(relxerr[te])))[round(nc-1-0.05*nc):nc-1]]
te=where(xm gt -1000 and ysig gt -1000 and ysig lt 4 and delta eq 1  and color ne cgcolor('darkgrey'),nc)
sry=te[((sort(relyerr[te])))[round(nc-1-0.05*nc):nc-1]]
terr=lonarr(n)&terr[[sry,srx]]=1

; Take the top 5% of x and y. Combine into array and plot both error bars. Only horizontal error bar if upper limit
Y_plot=ym_in & t3=where(delta eq 0) & Y_plot[t3]+=ysig[t3]


; upper limits
t2=where(xm gt -1000 and ysig gt -1000 and ysig lt 12 and delta eq 0 and terr eq 1,nc2);  and color gt -1,nc)
if nc2 gt 0 then $
for i=0,nc2-1 do begin
    	oplot,bhpred[t2[i]] +[-1,1]*xerr[t2[i]],[1,1]*Y_plot[t2[i]],color=long(color[t2[i]]),thick=0.7
    	;oplot,bhpred[t2[i]] *[1,1],[1,1]*ym_in[t2[i]]+[0,-0.5],color=long(color[t2[i]]),linestyle=1
    endfor

; error bars
t=where(xm gt -1000 and ysig gt -1000 and ysig lt 4 and delta eq 1 and terr eq 1,nc);  and color gt -1,nc)
for i=0,nc-1 do oploterror,bhpred[t[i]],ym_in[t[i]],xerr[t[i]],ysig[t[i]],/nohat,psym=3,color=long(color[t[i]]),thick=0.7


; DOTS
; upper limits
t2=where(xm gt -1000 and ysig gt -1000 and ysig lt 12 and delta eq 0,nc2);  and color gt -1,nc)
if nc2 gt 0 then $
	for i=0,nc2-1 do $
	    oplot,[bhpred[t2[i]]],[Y_plot[t2[i]]],psym=cgsymcat(11),color=long(color[t2[i]]),symsize=pointsize

t=where(xm gt -1000 and ysig gt -1000 and ysig lt 4 and delta eq 1,nc);  and color gt -1,nc)
for i=0,nc-1 do $
   oplot,[bhpred[t[i]]],[ym_in[t[i]]],psym=cgsymcat(16),color=long(color[t[i]]),symsize=pointsize


sigs=0.5 ; how far of the relation should names be allowed to be plotted

; first plot names below the relation.
c=alpha-epsilon*sigs
b=-1
a=beta[0]
dis=(beta[0]*bhpred-Y_plot+alpha-epsilon*sigs)/sqrt(beta[0]^2+1)
t=where(dis ge epsilon*0.1 and xm gt -1000 and ysig gt -1000 and ysig lt 12 ,nc )
if nc gt 6 then t=t[((sort(dis[t])))[nc-7:nc-1]]
if nc gt 0 then xyouts,bhpred[t]    ,Y_plot[t]    ,htname[t] ,align=-0.1,color=cgcolor('black'),charsize=!p.charsize/4
; Second  plot names above the relation.
c=alpha+epsilon*sigs
b=-1
a=beta[0]
dis=-(a*bhpred+b*Y_plot+c)/sqrt(a^2+b^2)
t=where(dis ge epsilon*0.1 and xm gt -1000 and ysig gt -1000 and ysig lt 12 and delta eq 1,nc)
if nc gt 6 then t=t[((sort(dis[t])))[nc-7:nc-1]]
if nc gt 0 then xyouts,bhpred[t]    ,Y_plot[t]    ,htname[t] ,align=1.1,color=cgcolor('black'),charsize=!p.charsize/4

textcolor=[cgcolor('green'),cgcolor('cyan'),cgcolor('red'),cgcolor('black'),cgcolor('darkgrey'),cgcolor('black')]
al_legend,['maser','gas','stars','reverberation','omitted','upper limit'],psym=[16,16,16,16,16,11],charsize=0.5 ,$
	 textcolor=textcolor,colors=textcolor,box=0 ,symsize=pointsize*1.25 

if !D.NAME eq 'PS' then device,/helvetica
sign='+'&if mean(post.beta[1]) lt 0 then sign = '-'
string=string(yunits+'=(',mean(post.alpha  ),'\pm',stddev(post.alpha  ),')+('             $
     ,mean(post.beta[0]),'\pm',stddev(post.beta[0]),')'+xunits[0]+sign+'('   $
	 ,abs(mean(post.beta[1])),'\pm',stddev(post.beta[1]),')'+xunits[1]+' (\epsilon=' $
	 ,mean(sqrt(post.sigsqr)),'\pm',stddev(sqrt(post.sigsqr)),')',$
		 format='(A,F5.1,A,F3.1,A,F4.2,A,F4.2,A,F4.2,A,F4.2,A,F4.2,A,F4.2,A)')
print,string
al_legend,textoidl(string),charsize=0.4,/right,box=0,/bottom
if keyword_set(file) then begin 
	device,/close
	spawn,'epstopdf '+file+' &'
	set_plot,'x'
endif

end

pro bh1Dlinmix,xm,xsig_in,ym_in,ysig,delta,s,htname,post=post,color,xrange=xrange,file=file,_EXTRA=_EXTRA,xunits=xunits,yunits=yunits
	if not keyword_set(xunits) then xunits=''
	if not keyword_set(yunits) then yunits=''
	
n=N_ELEMENTS(s)
linmix_err, xm[s], ym_in[s], post, Xsig=xsig_in[s], Ysig=ysig[s], delta=delta[s],/silent
alpha=mean(post.alpha)
beta=[mean(post.beta[0])]
epsilon=mean(sqrt(post.sigsqr))

bhpred=ALPHA + BETA ## xm 

if keyword_set(file) then begin
    set_plotstyle
	device,file=file
endif

t=where(xm gt -1000 and ysig gt -1000 and ysig lt 12 ,nc); and color gt -1,nc)
plot,xm[t],ym_in[t],/nodata,/yno,xrange=xrange,_EXTRA=_EXTRA ;$
t=where(xm gt -1000 and ysig gt -1000 and ysig lt 12 and delta eq 1,nc); and color gt -1,nc)

; plot the regression and the scatter.
oplot,!x.crange,ALPHA + BETA ##!x.crange                        ,color=cgcolor('black')
oplot,!x.crange,ALPHA + BETA ##!x.crange+epsilon,linestyle=2    ,color=cgcolor('grey')
oplot,!x.crange,ALPHA + BETA ##!x.crange-epsilon,linestyle=2    ,color=cgcolor('grey')
oplot,!x.crange,ALPHA + BETA ##!x.crange+epsilon*3,linestyle=1  ,color=cgcolor('grey')
oplot,!x.crange,ALPHA + BETA ##!x.crange-epsilon*3,linestyle=1  ,color=cgcolor('grey')


pointsize=0.4

; find out which error bars to plot
relxerr=xsig_in/abs(!x.crange[1]-!x.crange[0])
relyerr=ysig/abs(!y.crange[1]-!y.crange[0])
te=where(xm gt -1000 and ysig gt -1000 and  color ne cgcolor('darkgrey'),nc)
srx=te[((sort(relxerr[te])))[round(nc-1-0.05*nc):nc-1]]
te=where(xm gt -1000 and ysig gt -1000 and ysig lt 4 and delta eq 1  and color ne cgcolor('darkgrey'),nc)
sry=te[((sort(relyerr[te])))[round(nc-1-0.05*nc):nc-1]]
terr=lonarr(n)&terr[[sry,srx]]=1

; Take the top 5% of x and y. Combine into array and plot both error bars. Only horizontal error bar if upper limit.

; error bars of the upper limits 
t2=where(xm gt -1000 and ysig gt -1000 and ysig lt 12 and delta eq 0 and terr eq 1,nc2);  and color gt -1,nc)
if nc2 gt 0 then $
for i=0,nc2-1 do begin
    	oplot,xm[t2[i]] +[-1,1]*xsig_in[t2[i]],[1,1]*(ym_in+ysig)[t2[i]],color=long(color[t2[i]]),thick=0.7
    endfor

; error bars
t=where(xm gt -1000 and ysig gt -1000 and ysig lt 4 and delta eq 1 and terr eq 1,nc);  and color gt -1,nc)
for i=0,nc-1 do oploterror,xm[t[i]],ym_in[t[i]],xsig_in[t[i]],ysig[t[i]],/nohat,psym=3,color=long(color[t[i]]),thick=0.7


; DOTS
; upper limits 
t2=where(xm gt -1000 and ysig gt -1000 and ysig lt 12 and delta eq 0,nc2);  and color gt -1,nc)
if nc2 gt 0 then $
	for i=0,nc2-1 do $
	    oplot,[xm[t2[i]]],[(ym_in+ysig)[t2[i]]],psym=cgsymcat(11),color=long(color[t2[i]]),symsize=pointsize

t=where(xm gt -1000 and ysig gt -1000 and ysig lt 4 and delta eq 1,nc);  and color gt -1,nc)
for i=0,nc-1 do $
   oplot,[xm[t[i]]],[ym_in[t[i]]],psym=cgsymcat(16),color=long(color[t[i]]),symsize=pointsize


Y_plot=ym_in & t3=where(delta eq 0) & Y_plot[t3]+=ysig[t3]
dis=bhpred-y_plot
t=where(dis ge epsilon*0.1 and xm gt -1000 and ysig[s] gt -1000 and ysig lt 12  ,nc)
if nc gt 6 then t=t[((sort(dis[t])))[nc-7:nc-1]]
if nc gt 0 then xyouts,xm[t]    ,Y_plot[t]    ,htname[t] ,align=-0.1,color=cgcolor('black'),charsize=0.3
dis=y_plot-bhpred
t=where(dis ge epsilon*0.1 and xm gt -1000 and ysig gt -1000 and ysig lt 12 and delta eq 1 ,nc)
if nc gt 6 then t=t[((sort(dis[t])))[nc-7:nc-1]]
if nc gt 0 then xyouts,xm[t]    ,Y_plot[t]    ,htname[t] ,align=1.1,color=cgcolor('black'),charsize=0.3

textcolor=[cgcolor('green'),cgcolor('cyan'),cgcolor('red'),cgcolor('black'),cgcolor('darkgrey'),cgcolor('black')]
al_legend,['maser','gas','stars','reverberation','omitted','upper limit'],psym=[16,16,16,16,16,11],charsize=0.5 ,$
	 textcolor=textcolor,colors=textcolor,box=0 ,symsize=pointsize*1.25 

; Overplot the regression
string=string(yunits+'=(',mean(post.alpha  ),'\pm',stddev(post.alpha  ),')+('             $
	 ,(mean(post.beta)),'\pm',stddev(post.beta),')'+xunits+' (\epsilon=' $
	 ,mean(sqrt(post.sigsqr)),'\pm',stddev(sqrt(post.sigsqr)),')',$
		 format='(A,F5.1,A,F3.1,A,F4.2,A,F4.2,A,F4.2,A,F4.2,A)')
print,string
if !D.NAME eq 'PS' then device,/helvetica;,/narrow
al_legend,textoidl(string),charsize=0.4,/right,box=0,/bottom

if keyword_set(file) then begin 
	device,/close
	spawn,'epstopdf '+file+' &'
	set_plot,'x'
endif

end

