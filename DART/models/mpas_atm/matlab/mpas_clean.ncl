; This code is not protected by the DART copyright agreement.
; DART $Id$
 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

  r2d = 57.2957795             ; radians to degrees

  wks = gsn_open_wks("pdf","cells")
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

; TJH need the colormap to use in matlab
;  bob = gsn_retrieve_colormap(wks)
;  print(bob)

  f = addfile("output.nc","r")

  nEdgesOnCell   = f->nEdgesOnCell(:)
  verticesOnCell = f->verticesOnCell(:,:)
  lonCell        = f->lonCell(:) * r2d
  latCell        = f->latCell(:) * r2d
  lonVertex      = f->lonVertex(:) * r2d
  latVertex      = f->latVertex(:) * r2d

  res                      = True
  res@gsnPaperOrientation  = "portrait"

  res@sfXArray             = lonCell
  res@sfYArray             = latCell

  res@cnFillOn             = True
  res@cnFillMode           = "RasterFill"
  res@cnLinesOn            = False
  res@cnLineLabelsOn       = False
  res@cnInfoLabelOn        = False

  res@lbLabelAutoStride    = True
  res@lbBoxLinesOn         = False

  res@mpProjection      = "CylindricalEquidistant"
;  res@mpProjection      = "Orthographic"
  res@mpDataBaseVersion = "MediumRes"
  res@mpCenterLatF      = 0.
  res@mpCenterLonF      = 0.
  res@mpGridAndLimbOn   = False
  res@mpOutlineOn       = False
  res@mpFillOn          = False
  res@mpPerimOn         = False
  res@gsnFrame          = False

  ;
  ; The purpose of this section is simply to set up a graphic ('map')
  ;    that uses the projection specified above, and over which we
  ;    can draw polygons
  ;
  h      = f->areaCell(:)
  sizes  = dimsizes(h)
  nCells = sizes(0)
  xpoly  = new((/6/), "double")
  ypoly  = new((/6/), "double")
  res@cnConstFLabelOn = False
  res@lbLabelBarOn = False
  map = gsn_csm_contour_map(wks,h,res)

  ; This is the timestep index 
  t = stringtointeger(getenv("T"))
  t = 0
  copyindex = 0

  ;
  ; Set the field to be plotted here
  ;
  pres = True
  ; datmat   = f->uReconstructZonal(t,:,0)      for files without copies
  datmat   = f->uReconstructZonal(t,copyindex,:,0)
  minfld = min(datmat)
  maxfld = max(datmat)
  fldrange = maxfld - minfld
  do iCell=0,nCells-1
  do i=0,nEdgesOnCell(iCell)-1
     xpoly(i) = lonVertex(verticesOnCell(iCell,i)-1)
     ypoly(i) = latVertex(verticesOnCell(iCell,i)-1)
     if (i .gt. 0) then
        if (abs(xpoly(i) - xpoly(0)) .gt. 180.0) then
           if (xpoly(i) .gt. xpoly(0)) then
              xpoly(i) = xpoly(i) - 360.0
           else
              xpoly(i) = xpoly(i) + 360.0
           end if
        end if
     end if
  end do
  pres@gsFillColor = doubletointeger(198*(datmat(iCell) - minfld)/fldrange+2)
  gsn_polygon(wks,map,xpoly(0:nEdgesOnCell(iCell)-1),ypoly(0:nEdgesOnCell(iCell)-1),pres);
  end do


  ;
  ; Draw label bar
  ;

  xcb = new((/4/), "float")
  ycb = new((/4/), "float")

  tres = True
  tres@txAngleF = 90.0
  tres@txFontHeightF = 0.015
  do i=2,200
     xcb(0) = 0.125 + i*0.75/198
     ycb(0) = 0.11

     xcb(1) = 0.125 + (i+1)*0.75/198
     ycb(1) = 0.11

     xcb(2) = 0.125 + (i+1)*0.75/198
     ycb(2) = 0.16

     xcb(3) = 0.125 + i*0.75/198
     ycb(3) = 0.16

     tres@gsFillColor = i

     gsn_polygon_ndc(wks,xcb,ycb,tres);

     j = (i-2) % 20
     if ((j .eq. 0) .or. (i .eq. 200)) then
        ff = minfld + int2flt(i-2) * fldrange / 198.0
        label = sprintf("%5.3g", ff)
        gsn_text_ndc(wks, label, xcb(0), 0.060, tres)
     end if

  end do

  mres = True
  mres@mpCenterLatF      = 0.
  mres@mpCenterLonF      = 0.
  mres@mpGridAndLimbOn   = False

; turn mpOutlineOn False for no coastlines, True for coastlines
  mres@mpOutlineOn       = False
  mres@mpFillOn          = False
  mres@mpPerimOn         = False
  mres@gsnFrame          = False
  mapo = gsn_csm_map(wks,mres)

  frame(wks)

end

; <next few lines under version control, do not edit>
; $URL$
; $Revision$
; $Date$
 
