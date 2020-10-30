; Program started June 2015
; Explore dotwnloaded SDSS data for pea sample
; last edited 6/29/2017

;Part 0 Reading in & Making calculations about the original 251/80 sample
;--------------------------------------------------------
;Get file with Pea data from original pea search 2009, DR7:
   ; file saved finalpea_ccardamone.fit, in cat
   readcol,'cat/2009_PEAFLAGS_ccardamone.csv',objid,ra,raerr,dec,decerr,type,psfmag_u,psfmagErr_u,psfmag_g,psfmagErr_g,psfmag_r,psfmagErr_r,psfmag_i,psfmagErr_i,psfmag_z,psfmagErr_z,extinction_u,extinction_g,extinction_r,extinction_i,extinction_z,z,zerr,zConf,SpecClass,objTypename,run,rerun,camcol,field,obj,mjd,plate,fiberID,sn_0,sn_1,sn_2,velDisp,velDispErr,petroRad_r,ew,format='(A,D,D,D,D, F, F,F,F,F,F,F,F,F,F,F, F,F,F,F,F, F,F,F,I,A, I,I,I,I,I,L,I,I, F,F,F,F,F,F,F)'
   umag=psfmag_u-extinction_u
   gmag=psfmag_g-extinction_g
   rmag=psfmag_r-extinction_r
   imag=psfmag_i-extinction_i
   zmag=psfmag_z-extinction_z
   
  ;& 80 objects confirmed re spectra as SB
   readcol,'cat/2009_SB_data.txt',SDSSID,arrloc,redshift,ebmv,ha_sfr,ha_sfr_cor,o2_sfr,nuv_sfr,fuv_sfr,o3_sig,sig_stars,format='(A,I,F,F,F,F,F,F,F,F,F)'
;   ;Luminous Red Galaxy Samples
;    ; BOSS LRG LOWZ sample from Nikhil
;                                ; (https://data.sdss.org/sas/dr12/boss/lss/, https://data.sdss.org/datamodel/files/BOSS_LSS_REDUX/
;;   ;created from  galaxy_DR12v5_LOWZ_North.fits and    galaxy_DR12v5_LOWZ_South.fits
;     LRG=mrdfits('cat/LRGzcat.fits',1)  


; AFTER running colorselect cats to combine catalogs from SDSS SQL
; querys with full color criteria and to match to portsmouth catalogs
; spectra by position....

   ;   ;Note - * Portsmouth Group Galaxy Fitting w/ GANDALF (pop models for continuum from Maraston 2011, THomas 2011Classification from Kauffmann et al. (2003), Kewley et al. (2001) and Schawinski et al. (2007) for objects with the emission lines H&amp;beta;, [OIII], H&amp;alpha;, [NII] available with A;/N &gt; 1.5. Possible values: "BLANK" (if emission lines not available), "Star Forming", "Seyfert", "LINER", "Seyfert/LINER", "Composite"
;;http://www.sdss.org/dr12/spectro/galaxy_portsmouth/
   ;;; NOTE these files were created by position matching to the
   ;;; portsmouth catalog
       dr12peaport=mrdfits('output/dr12peaport.fits',1)
       dr14peaport=mrdfits('output/dr14peaport.fits',1)
   ;;;; thesea re the full color search without any matching to
   ;;;; portsmouth data
       dr12pea=mrdfits('cat/peas2017dr12fullcolor_noport.fits',1)
       dr14pea=mrdfits('cat/peas2017dr14fullcolor_noport.fits',1)
   ;;;;run dr12peascatalog from here with the full color, through
   ;;;;CasJobs to hunt dr 14 nieghbors with 10 arcseconds....
       dr14mdr12=mrdfits('output/dr14matchdr12peas.fits',1)

;;;;Note in 2009 sample, I rejected ~5 objects, unfittable spectra,
;;;;skyline over OIII and stange bpt location

       norayesarr=fltarr(N_ELEMENTS(dr14pea.objid))-1
       norayesarrport=fltarr(n_Elements(dr14peaport.objid))-1
       readcol,'output/norayes.cat',noraobjid,norara,noradec,noraz,format='(A,D,D,F)'
       for i=0,n_elements(noraobjid)-1 do begin
          s=where(dr14pea.objid eq noraobjid(i))
;          print, noraobjid(i),dr14pea(s).objid,norara(i),dr14pea(s).ra
          norayesarr(s)=1
          sport=where(dr14peaport.objid eq noraobjid(i))
          if sport ne -1 then norayesarrport(sport)=1
;          print,i, sport,noraobjid(i),format='(I,3x,A,3x,A)'
       endfor
;       
       
       ;;; Question - why some different when psoiton match vs. catlog
       ;;;            id match in dr 14
dr14CasJobport=mrdfits('cat/peas2017dr14.fits',1)

;So which aren't retained in DR14 but were in dr 12 and WHY did we loose some?
  for i=0,n_elements(dr12pea)-1 do begin
     gcirc,1,dr12pea[i].ra/15,dr12pea[i].dec,dr14pea.ra/15,dr14pea.dec,d2
     dmatch=min(d2,m)  ; assigns min dist to dmatch, its location to m
                       ;distance  is in arcsec,
    if i eq 0 then dr12dr14match=m else  dr12dr14match=[dr12dr14match,m]
    if i eq 0 then dr12dr14dist=dmatch else dr12dr14dist=[dr12dr14dist,dmatch]    
 endfor

  tgood=dr12dr14match
  tgood[where(dr12dr14dist gt 1)]=-1

  t12=where(dr12dr14dist lt 1)
  t14=dr12dr14match(where(dr12dr14dist lt 1))

  print, 'peas with same RA/Dec in DR12 v. DR14 but different redshift'
  
  tdiff=where( (abs(dr12pea.z-dr14pea[dr12dr14match].z) gt .001) and (dr12dr14dist lt 1))
  for i = 0,n_elements(tdiff)-1 do begin
   print, dr12pea(tdiff(i)).objid,dr12pea(tdiff(i)).ra,dr12pea(tdiff(i)).dec
     endfor
  
;  print, 'DR12 peas w/out a dr14 pea match', 
;;;Note these distances are either ~0 (or ~.2 for rounding in dr7 data
;;;above), or very large >> 100 arcseconds

t=where(dr12dr14dist gt 1)
 print, 'PEAS not in DR14 but in dr12 catalog'
  for i=0,n_elements(t)-1 do print, dr12pea[t[i]].objid,dr12pea[t[i]].ra,dr12pea[t[i]].dec

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Investigate sample - look at color-color plots
  qso=mrdfits('../peas/cat/qsocomp.fit',1)
  gal=mrdfits('../peas/cat/galcompewflag_ccardamone.fit',1)
   match=mrdfits('../peas/cat/matcharr.fits',1)

;note if tflag=yes, then where tgood ne -1 lines connect cat 1 and
;cat2 in color plot
    l1='Cardamone et al. 2009'   ;
    l2='Nora Confirmed Peas'  ;
    l3='All Peas  DR14'
    lgal='C09 Comparison Galaxy Sample'
    lqso='C09 Comparsion QSO Sample'
;;;;; Color-color plot codes
psfile='plots/dr1214color_gr_ri.ps'                                         
    xr=[-1.8,.75]
    yr=[-.5,2.5]
    xtitle='r-i'
    ytitle='g-r'
    xval=rmag-imag
    yval=gmag-rmag

    xval1=(dr12pea.psfmag_r-dr12pea.extinction_r)-(dr12pea.psfmag_i-dr12pea.extinction_i)
    yval1=(dr12pea.psfmag_g-dr12pea.extinction_g)-(dr12pea.psfmag_r-dr12pea.extinction_r)
   tgood=indgen(n_elements(Xval1)) 
   t=(where(norayesarr eq 1))
   xval2=(dr14pea[t].psfmag_r-dr14pea[t].extinction_r)-(dr14pea[t].psfmag_i-dr14pea[t].extinction_i)
   yval2=(dr14pea[t].psfmag_g-dr14pea[t].extinction_g)-(dr14pea[t].psfmag_r-dr14pea[t].extinction_r)
    xval3=(dr14pea[t].psfmag_r-dr14pea[t].extinction_r)-(dr14pea[t].psfmag_i-dr14pea[t].extinction_i)
    yval3=(dr14pea[t].psfmag_g-dr14pea[t].extinction_g)-(dr14pea[t].psfmag_r-dr14pea[t].extinction_r)
    galxval=gal(match.gal).psfmag_r-gal(match.gal).extinction_r-(gal(match.gal).psfmag_i-gal(match.gal).extinction_i)
    galyval=gal(match.gal).psfmag_g-gal(match.gal).extinction_g-(gal(match.gal).psfmag_r-gal(match.gal).extinction_r)
    qsoxval=qso.psfmag_r-qso.extinction_r-(qso.psfmag_i-qso.extinction_i)
    qsoyval=qso.psfmag_g-qso.extinction_g-(qso.psfmag_r-qso.extinction_r)
;;;tgood=1
    collab='gr_ri'
   tflag='yes'
    mkpeacolplot,psfile,xtitle,ytitle,xr,yr,xval,yval,xval1,yval1,xval2,yval2,xval3,yval3,galxval,galyval,qsoxval,qsoyval,collab,l1,l2,l3,lgal,lqso,tgood,tflag

psfile='plots/dr12dr14color_ur_rz.ps'                                         
    xr=[-1.6,1.2]
    yr=[-.5,4]
    xtitle='r-z'
    ytitle='u-r'
    xval=rmag-zmag
    yval=umag-rmag
    xval1=(dr12pea.psfmag_r-dr12pea.extinction_r)-(dr12pea.psfmag_z-dr12pea.extinction_z)
    yval1=(dr12pea.psfmag_u-dr12pea.extinction_u)-(dr12pea.psfmag_r-dr12pea.extinction_r)
tgood=indgen(n_elements(Xval1))
    xval2=(dr14mdr12.psfmag_r-dr14mdr12.extinction_r)-(dr14mdr12.psfmag_z-dr14mdr12.extinction_z)
    yval2=(dr14mdr12.psfmag_u-dr14mdr12.extinction_u)-(dr14mdr12.psfmag_r-dr14mdr12.extinction_r)
    xval3=(dr14pea.psfmag_r-dr14pea.extinction_r)-(dr14pea.psfmag_z-dr14pea.extinction_z)
    yval3=(dr14pea.psfmag_u-dr14pea.extinction_u)-(dr14pea.psfmag_r-dr14pea.extinction_r)

    galxval=gal(match.gal).psfmag_r-gal(match.gal).extinction_r-(gal(match.gal).psfmag_z-gal(match.gal).extinction_z)
    galyval=gal(match.gal).psfmag_u-gal(match.gal).extinction_u-(gal(match.gal).psfmag_r-gal(match.gal).extinction_r)
    qsoxval=qso.psfmag_r-qso.extinction_r-(qso.psfmag_z-qso.extinction_z)
    qsoyval=qso.psfmag_u-qso.extinction_u-(qso.psfmag_r-qso.extinction_r)


    collab='ur_rz'
 mkpeacolplot,psfile,xtitle,ytitle,xr,yr,xval,yval,xval1,yval1,xval2,yval2,xval3,yval3,galxval,galyval,qsoxval,qsoyval,collab,l1,l2,l3,lgal,lqso,tgood,tflag

 


    STOP
 end


;note if tflag=yes, then where tgood ne -1 lines connect cat 1 and
;cat2 in color plot
    l1='DR12'   ;
    l2='DR14 colors of dr12 selection'  ;
    l3='DR14'
    lgal='C09 Comparison Galaxy Sample'
    lqso='C09 Comparsion QSO Sample'
;;;;; Color-color plot codes
psfile='plots/dr1214color_gr_ri.ps'                                         
    xr=[-1.8,.75]
    yr=[-.5,2.5]
    xtitle='r-i'
    ytitle='g-r'
    xval=rmag-imag
    yval=gmag-rmag

    xval1=(dr12pea.psfmag_r-dr12pea.extinction_r)-(dr12pea.psfmag_i-dr12pea.extinction_i)
    yval1=(dr12pea.psfmag_g-dr12pea.extinction_g)-(dr12pea.psfmag_r-dr12pea.extinction_r)
tgood=indgen(n_elements(Xval1))
    xval2=(dr14mdr12.psfmag_r-dr14mdr12.extinction_r)-(dr14mdr12.psfmag_i-dr14mdr12.extinction_i)
    yval2=(dr14mdr12.psfmag_g-dr14mdr12.extinction_g)-(dr14mdr12.psfmag_r-dr14mdr12.extinction_r)
    xval3=(dr14pea.psfmag_r-dr14pea.extinction_r)-(dr14pea.psfmag_i-dr14pea.extinction_i)
    yval3=(dr14pea.psfmag_g-dr14pea.extinction_g)-(dr14pea.psfmag_r-dr14pea.extinction_r)

    galxval=gal(match.gal).psfmag_r-gal(match.gal).extinction_r-(gal(match.gal).psfmag_i-gal(match.gal).extinction_i)
    galyval=gal(match.gal).psfmag_g-gal(match.gal).extinction_g-(gal(match.gal).psfmag_r-gal(match.gal).extinction_r)
    qsoxval=qso.psfmag_r-qso.extinction_r-(qso.psfmag_i-qso.extinction_i)
    qsoyval=qso.psfmag_g-qso.extinction_g-(qso.psfmag_r-qso.extinction_r)
;;;tgood=1
    collab='gr_ri'
   tflag='yes'
    mkpeacolplot,psfile,xtitle,ytitle,xr,yr,xval,yval,xval1,yval1,xval2,yval2,xval3,yval3,galxval,galyval,qsoxval,qsoyval,collab,l1,l2,l3,lgal,lqso,tgood,tflag

psfile='plots/dr12dr14color_ur_rz.ps'                                         
    xr=[-1.6,1.2]
    yr=[-.5,4]
    xtitle='r-z'
    ytitle='u-r'
    xval=rmag-zmag
    yval=umag-rmag
    xval1=(dr12pea.psfmag_r-dr12pea.extinction_r)-(dr12pea.psfmag_z-dr12pea.extinction_z)
    yval1=(dr12pea.psfmag_u-dr12pea.extinction_u)-(dr12pea.psfmag_r-dr12pea.extinction_r)
tgood=indgen(n_elements(Xval1))
    xval2=(dr14mdr12.psfmag_r-dr14mdr12.extinction_r)-(dr14mdr12.psfmag_z-dr14mdr12.extinction_z)
    yval2=(dr14mdr12.psfmag_u-dr14mdr12.extinction_u)-(dr14mdr12.psfmag_r-dr14mdr12.extinction_r)
    xval3=(dr14pea.psfmag_r-dr14pea.extinction_r)-(dr14pea.psfmag_z-dr14pea.extinction_z)
    yval3=(dr14pea.psfmag_u-dr14pea.extinction_u)-(dr14pea.psfmag_r-dr14pea.extinction_r)

    galxval=gal(match.gal).psfmag_r-gal(match.gal).extinction_r-(gal(match.gal).psfmag_z-gal(match.gal).extinction_z)
    galyval=gal(match.gal).psfmag_u-gal(match.gal).extinction_u-(gal(match.gal).psfmag_r-gal(match.gal).extinction_r)
    qsoxval=qso.psfmag_r-qso.extinction_r-(qso.psfmag_z-qso.extinction_z)
    qsoyval=qso.psfmag_u-qso.extinction_u-(qso.psfmag_r-qso.extinction_r)


    collab='ur_rz'
 mkpeacolplot,psfile,xtitle,ytitle,xr,yr,xval,yval,xval1,yval1,xval2,yval2,xval3,yval3,galxval,galyval,qsoxval,qsoyval,collab,l1,l2,l3,lgal,lqso,tgood,tflag

 


    STOP
 end



;;; Make Plots of samples
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Create Plots . . .
;   ;Specify the red component of each color:
;    RED = [0, 1, 1, 0, 0, 1]
;    ;Specify the green component of each color:
;    GREEN = [0, 1, 0, 1, 0, 1]
;    ;Specify the blue component of each color:
;    BLUE = [0, 1, 0, 0, 1, 0]
;     ;Load the first six elements of the color table:
;    TVLCT, 255 * RED, 255 * GREEN, 255 * BLUE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;    colorarr=[fsc_color('Cyan',!D.Table_Size-2),fsc_color('Dark Green',!D.Table_Size-3),$
;              fsc_color('Spring Green',!D.Table_Size-4),fsc_color('Steel Blue',!D.Table_Size-5),$
;              fsc_color('Dark Orchid',!D.Table_Size-6),fsc_color('Orchid',!D.Table_Size-7),$
;              fsc_color('Red',!D.Table_Size-8),fsc_color('Navy',!D.Table_Size-9)]
;    ;Print, FSC_Color(/Names), Format='(6A15)'
;
;     gal=mrdfits('../peas/cat/galcompewflag_ccardamone.fit',1)
;     match=mrdfits('../peas/cat/matcharr.fits',1)
;
;     
;
;    
;;    
;;    l1='Original 80 Peas'        ; original 80 solid sample
;;    l2='DR12 Original Peas Recovered'   ; all star forming objects
;;    l3='DR12 color-sel SF Peas 479 obj'  ;all starforming + good redshifts
;;    l4='Luminous Red Galaxies (200k+)' ;LRGs, not this only counts those with 0<z<0.4
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Create Plots . . .
;; Plot redshift distributions 
;    psfile='plots/pea_zdist.ps'
;    binsize=0.01
;    peazhist=histogram(z,locations=lpeaz,min=0,max=.4,binsize=binsize)  ;E251 sample, all
;    sfpeazhist=histogram(z(arrloc),locations=lsfpeaz,min=0,max=.4,binsize=binsize)
;    newpeahist=histogram(peaspec[sfspec].z,locations=lnewpeaz,min=0,max=.4,binsize=binsize)
;    newpeagood1hist=histogram(peaspec[sfspecs].z,locations=lnewgood1peaz,min=0,max=.4,binsize=binsize)
;    lrgzhist=histogram(lrg.z,locations=llrg,min=0,max=.4,binsize=binsize) 
;    xr=[0,0.4]
;    yr=[0,0.12]
;    xtitle='redshift'
;    ytitle='Normalized Histogram'
;    xval=lpeaz                  ; left hand side of bin (note plot center when using psym=10)
;    yval=sfpeazhist/total(sfpeazhist)
;    yval2=newpeahist/total(newpeahist)
;    yval3=newpeagood1hist/total(newpeagood1hist)
;    yval4=lrgzhist/total(lrgzhist)
;    lx=findgen(10)*.001+0.02
;    ly=fltarr(10)+.11
;    ly2=fltarr(10)+.105
;    ly3=fltarr(10)+.10
;    ly4=fltarr(10)+.095
;    ylog='no'
;
;mkpeaplot,psfile,xtitle,ytitle,xr,yr,xval,yval,yval2,yval3,yval4,lx,ly,ly2,ly3,ly4,l1,l2,l3,l4,binsize,ylog
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
;; Make histograms & Plot of the Equivalent Width of OIII / o3
;
;  psfile='plots/o3ewhist.ps'
;  binsize=50
;    peaew=histogram(ew[arrloc],min=-10,max=2500,locations=lpeaew,binsize=binsize)
;    newpeaew=histogram(peaspec[sfspec].EW_OIII_5006,min=-10,max=2500,locations=lpea,binsize=binsize)
;    newpeagood1ew=histogram(peaspec[sfspecs].EW_OIII_5006,min=-10,max=2500,locations=lpea,binsize=binsize)
;    galew=histogram(gal(match.gal).ew,min=-10,max=2500,locations=lgal,binsize=binsize)
;    xr=[-1,1800]
;    yr=[1,1000]
;    xtitle='[O III] Equivalent Width [A]'
;    ytitle='Histogram'
;    xval=lpeaew          ; left hand side of bin (note plot center when using psym=10)
;    yval=peaew
;    yval2=newpeaew
;    yval3=newpeagood1ew
;    yval4=galew
;    lx=findgen(20)+200
;    ly=fltarr(20)+700
;    ly2=fltarr(20)+450
;    ly3=fltarr(20)+325
;    ly4=fltarr(20)+250
;    l4='Galaxy Comparison Sample C09 (~12k)'
;    ylog='yes'
;    
;    mkpeaplot,psfile,xtitle,ytitle,xr,yr,xval,yval,yval2,yval3,yval4,lx,ly,ly2,ly3,ly4,l1,l2,l3,l4,binsize,ylog
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
;; Make histograms & Plot of size of the galaxies (Petrorad_R in kpc
;    
;   psfile='plots/petroradhist.ps'
;   binsize = 1 ;kpc
;  ;;uses z and gal.z to compute angle
;  petrodkpc=zang(1,z) ;distance in arcseconds of 1 KPC at redshift of pea
;  peahistkpc=histogram(petrorad_r(arrloc)/petrodkpc(arrloc),locations=lpeakpc,min=0,max=100,binsize=binsize)
;  newpetrodkpc=zang(1,peaspec[sfspec].z)
;  newpeahistkpc=histogram(peaspec[sfspecs].petrorad_r/petrodkpc[sfclean],locations=lpeakpc1,min=0,max=100,binsize=binsize)
;  good1newpeahistkpc=histogram(peaspec[sfspecs].petrorad_r/petrodkpc[good1],locations=lpeakpc2,min=0,max=100,binsize=binsize)
;  galpetrodkpc=zang(1,gal.z)
;  galhistkpc=histogram(gal(match.gal).petrorad_r/galpetrodkpc(match.gal),locations=lgalkpc,min=0,max=100,binsize=binsize)
;
;  
;;kstwo,peahist,galhist,D,prob
;;print, 'KS test results ["]: prob=',prob,'d=',d
;;kstwo,peahistkpc,galhistkpc,D,prob
;;print, 'KS test results [Kpc]: prob=',prob,'d=',d;;;;
;
;    xr=[0,50]
;    yr=[.9,1000]
;    xtitle='Petrorad_R in kpc'
;    ytitle='Histogram'
;    xval=lpeakpc          ; left hand side of bin (note plot center when using psym=10)
;    yval=peahistkpc
;    yval2=newpeahistkpc
;    yval3=good1newpeahistkpc
;    yval4=galhistkpc
;    lx=findgen(2)+30
;    ly=fltarr(2)+700
;    ly2=fltarr(2)+450
;    ly3=fltarr(2)+325
;    ly4=fltarr(2)+250
;    l4='Galaxy Comparison Sample C09 (~12k)'
;    ylog='yes'
;    
;    mkpeaplot,psfile,xtitle,ytitle,xr,yr,xval,yval,yval2,yval3,yval4,lx,ly,ly2,ly3,ly4,l1,l2,l3,l4,binsize,ylog



;;Investigate sample - look at color-color plots, OIII dists, etc.
;  qso=mrdfits('../peas/cat/qsocomp.fit',1)
;    l1='Original 250 Peas'        ; original 80 solid sample
;    l2='DR12 Spec nonSF'   ; all star forming objects
;    l3='DR14 Spec nonSF'  ;all starforming + good redshifts
;    lgal='C09 Comparison Galaxy Sample'
;    lqso='C09 Comparsion QSO Sample'
;;;;; Color-color plot codes
;psfile='plots/color_gr_ri.ps'                                         
;
;    xr=[-1.8,.75]
;    yr=[-.5,2.5]
;    xtitle='r-i'
;    ytitle='g-r'
;    xval=rmag-imag
;    yval=gmag-rmag
;    ;this was for the 80 BPT selected
;;    xval1=rmag(arrloc)-imag(arrloc)
;;    yval1=gmag(arrloc)-rmag(arrloc)
;;    xval2=(peaspec[peaspecorig].psfmag_r-peaspec[peaspecorig].extinction_r)-(peaspec[peaspecorig].psfmag_i-peaspec[peaspecorig].extinction_i)
;;    yval2=(peaspec[peaspecorig].psfmag_g-peaspec[peaspecorig].extinction_g)-(peaspec[peaspecorig].psfmag_r-peaspec[peaspecorig].extinction_r)
;;    xval3=(peaspec[sfspecs].psfmag_r-peaspec[sfspecs].extinction_r)-(peaspec[sfspecs].psfmag_i-peaspec[sfspecs].extinction_i)
;;    yval3=(peaspec[sfspecs].psfmag_g-peaspec[sfspecs].extinction_g)-(peaspec[sfspecs].psfmag_r-peaspec[sfspecs].extinction_r)
;
;
;    xval1=rmag-imag
;    yval1=gmag-rmag
;
;    xval2=(dr12[seldr12].psfmag_r-dr12[seldr12].extinction_r)-(dr12[seldr12].psfmag_i-dr12[seldr12].extinction_i)
;    yval2=(dr12[seldr12].psfmag_g-dr12[seldr12].extinction_g)-(dr12[seldr12].psfmag_r-dr12[seldr12].extinction_r)
;    xval3=(dr14[seldr14].psfmag_r-dr14[seldr14].extinction_r)-(dr14[seldr14].psfmag_i-dr14[seldr14].extinction_i)
;    yval3=(dr14[seldr14].psfmag_g-dr14[seldr14].extinction_g)-(dr14[seldr14].psfmag_r-dr14[seldr14].extinction_r)
;
;    galxval=gal(match.gal).psfmag_r-gal(match.gal).extinction_r-(gal(match.gal).psfmag_i-gal(match.gal).extinction_i)
;    galyval=gal(match.gal).psfmag_g-gal(match.gal).extinction_g-(gal(match.gal).psfmag_r-gal(match.gal).extinction_r)
;    qsoxval=qso.psfmag_r-qso.extinction_r-(qso.psfmag_i-qso.extinction_i)
;    qsoyval=qso.psfmag_g-qso.extinction_g-(qso.psfmag_r-qso.extinction_r)
;tgood=1
;    collab='gr_ri'
;    mkpeacolplot,psfile,xtitle,ytitle,xr,yr,xval,yval,xval1,yval1,xval2,yval2,xval3,yval3,galxval,galyval,qsoxval,qsoyval,collab,l1,l2,l3,lgal,lqso,tgood
;    
;psfile='plots/color_ur_rz.ps'                                         
;
;    xr=[-1.6,1.2]
;    yr=[-.5,4]
;    xtitle='r-z'
;    ytitle='u-r'
;    xval=rmag-zmag
;    yval=umag-rmag
;    xval1=rmag(arrloc)-zmag(arrloc)
;    yval1=umag(arrloc)-rmag(arrloc)
;
;;    xval2=(peaspec[peaspecorig].psfmag_r-peaspec[peaspecorig].extinction_r)-(peaspec[peaspecorig].psfmag_z-peaspec[peaspecorig].extinction_z)
;;    yval2=(peaspec[peaspecorig].psfmag_u-peaspec[peaspecorig].extinction_u)-(peaspec[peaspecorig].psfmag_r-peaspec[peaspecorig].extinction_r)
;;    xval3=(peaspec[sfspecs].psfmag_r-peaspec[sfspecs].extinction_r)-(peaspec[sfspecs].psfmag_z-peaspec[sfspecs].extinction_z)
;;    yval3=(peaspec[sfspecs].psfmag_u-peaspec[sfspecs].extinction_u)-(peaspec[sfspecs].psfmag_r-peaspec[sfspecs].extinction_r)
;
;    xval2=(dr12[seldr12].psfmag_r-dr12[seldr12].extinction_r)-(dr12[seldr12].psfmag_z-dr12[seldr12].extinction_z)
;    yval2=(dr12[seldr12].psfmag_u-dr12[seldr12].extinction_u)-(dr12[seldr12].psfmag_r-dr12[seldr12].extinction_r)
;    xval3=(dr14[seldr14].psfmag_r-dr14[seldr14].extinction_r)-(dr14[seldr14].psfmag_z-dr14[seldr14].extinction_z)
;    yval3=(dr14[seldr14].psfmag_u-dr14[seldr14].extinction_u)-(dr14[seldr14].psfmag_r-dr14[seldr14].extinction_r)
;
;    galxval=gal(match.gal).psfmag_r-gal(match.gal).extinction_r-(gal(match.gal).psfmag_z-gal(match.gal).extinction_z)
;    galyval=gal(match.gal).psfmag_u-gal(match.gal).extinction_u-(gal(match.gal).psfmag_r-gal(match.gal).extinction_r)
;    qsoxval=qso.psfmag_r-qso.extinction_r-(qso.psfmag_z-qso.extinction_z)
;    qsoyval=qso.psfmag_u-qso.extinction_u-(qso.psfmag_r-qso.extinction_r)
;
;
;    collab='ur_rz'
;    mkpeacolplot,psfile,xtitle,ytitle,xr,yr,xval,yval,xval1,yval1,xval2,yval2,xval3,yval3,galxval,galyval,qsoxval,qsoyval,collab,l1,l2,l3,lgal,lqso,tgood
;
;
;;plot,rmag-zmag,umag-rmag,psym=3,xtitle='r-z',ytitle='u-r',xrange=[-1.6,1.2],$
;;   yrange=[-.5,4],/ystyle,/xstyle,/nodata,charsize=1.5,charthick=2;
;
;;;; OUtput to a file the list of peas, SDSSid,ra,dec,z,zerr,photflag, magnitude
;output_filename='output/table_DR12PeasRADEC.txt'
; print, 'outputfile='+output_filename
; openw,lun,/GET_LUN,output_filename
; printf,lun,'OBJID       RA       DEC'
; 
;    for i=0,n_elements(sfspec)-1 do begin
;       printf,lun,peaspec[sfspec[i]].objid,peaspec[sfspec[i]].ra,peaspec[sfspec[i]].dec, $
;              format='(A18,F12.6, F12.6)'
;endfor
;    free_lun,lun
;
;sellrg=where(lrg.z gt 0.112 and lrg.z lt 0.36)
;
;    
;    output_filename='output/table_exploreDR12_photflag.txt'
; print, 'outputfile='+output_filename
; openw,lun,/GET_LUN,output_filename
; printf,lun,'OBJID       RA       DEC'
; 
;    for i=0,n_elements(sfnot)-1 do begin
;       printf,lun,peaspec[sfnot[i]].objid,peaspec[sfnot[i]].ra,peaspec[sfnot[i]].dec, $
;              format='(A18,F12.6, F12.6)'
;endfor
;    free_lun,lun
;
; output_filename='output/table_DR13_colorz.txt'
; print, 'outputfile='+output_filename
; openw,lun,/GET_LUN,output_filename
; printf,lun,'     DR12_OBJID       ra  dec   z'
; 
; 
; 
;    for i=0,n_elements(seldr13col)-1 do begin
;       printf,lun,dr13p[seldr13col(i)].OBJID,dr13p[seldr13col(i)].ra,dr13p[seldr13col(i)].dec,dr13p[seldr13col(i)].z,$
;              format='(A18,F12.6, F12.6, F8.4)'
;endfor
;    free_lun,lun
;
;
;    output_filename='output/table_DR12colselnonSeyfert.txt'
; print, 'outputfile='+output_filename
; openw,lun,/GET_LUN,output_filename
; printf,lun,'          ID          ra    dec   z'
; 
;    for i=0,n_elements(selDR12nonSeyfert)-1 do begin
;       printf,lun,peaspec[selDR12nonSeyfert(i)].OBJID,peaspec[selDR12nonSeyfert(i)].ra,peaspec[selDR12nonSeyfert(i)].dec,peaspec[selDR12nonSeyfert(i)].z,$
;              format='(A18,F12.6, F12.6, F8.4)'
;endfor
;    free_lun,lun
;
;
;    
; END
;
;
;;;;;;;;;;;;;;;;
; 
;
;output_filename='output/table_NewPeaSFcleanphot.txt'
; print, 'outputfile='+output_filename
; openw,lun,/GET_LUN,output_filename
; printf,lun,'     DR12_OBJID       ra  dec   z'
; 
;    for i=0,n_elements(sfclean)-1 do begin
;       printf,lun,peaspec[sfclean(i)].OBJID,peaspec[sfclean(i)].ra,peaspec[sfclean(i)].dec,peaspec[sfclean(i)].z,$
;              format='(A18,F12.6, F12.6, F8.4)'
;endfor
;    free_lun,lun
;
;output_filename='output/table_NewPeaSFcleanphot.txt'
; print, 'outputfile='+output_filename
; openw,lun,/GET_LUN,output_filename
; printf,lun,'     DR12_OBJID       ra  dec   z'
; 
;    for i=0,n_elements(sfclean)-1 do begin
;       printf,lun,peaspec[sfclean(i)].OBJID,peaspec[sfclean(i)].ra,peaspec[sfclean(i)].dec,peaspec[sfclean(i)].z,$
;              format='(A18,F12.6, F12.6, F8.4)'
;endfor
;    free_lun,lun
;
;
;end
;
