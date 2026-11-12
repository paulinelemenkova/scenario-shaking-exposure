#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
NAME=fig_site_effect
mkdir -p work_site && cd work_site
CITIES=(kabul herat mazar kunduz jalalabad kandahar)
fetch() { if [[ $1 == http* ]]; then curl -s -f -o $2 $1/$2 || true; else cp $1/$2 . 2>/dev/null || true; fi; }
for c in "${CITIES[@]}"; do
  fetch $SITE ${c}_ratio.tif
  fetch $PREP ${c}_built.tif; fetch $PREP ${c}_alluvium.gmt; fetch $PREP ${c}_fan.gmt
done
fetch $SITE site_effect.csv

declare -A LON LAT TITLE
LON=([kabul]=69.181 [herat]=62.200 [mazar]=67.110 [kunduz]=68.860 [jalalabad]=70.450 [kandahar]=65.710)
LAT=([kabul]=34.519 [herat]=34.348 [mazar]=36.709 [kunduz]=36.728 [jalalabad]=34.430 [kandahar]=31.620)
TITLE=([kabul]=Kabul [herat]=Herat [mazar]=Mazar-i-Sharif [kunduz]=Kunduz [jalalabad]=Jalalabad [kandahar]=Kandahar)
LET=(a b c d e f)
PW=5.05; GAP=1.35; HX=0.30; HY=0.24

cat > ratio.cpt <<'CPT'
0.80 33/102/172  0.90 33/102/172  ;0.80-0.90
0.90 146/197/222 0.97 146/197/222 ;0.90-0.97
0.97 247/247/247 1.03 247/247/247 ;0.97-1.03
1.03 253/219/199 1.10 253/219/199 ;1.03-1.10
1.10 244/165/130 1.20 244/165/130 ;1.10-1.20
1.20 214/96/77   1.35 214/96/77   ;1.20-1.35
1.35 178/24/43   1.60 178/24/43   ;1.35-1.60
B 33/102/172
F 178/24/43
N 255/255/255
CPT
gmt makecpt -Cgray -T-1/1 -G0.55/1 > shade.cpt

for c in "${CITIES[@]}"; do
  lo=${LON[$c]}; la=${LAT[$c]}
  R=$(awk -v x=$lo -v y=$la -v hx=$HX -v hy=$HY 'BEGIN{printf "%.4f/%.4f/%.4f/%.4f", x-hx, x+hx, y-hy, y+hy}')
  echo $R > R_$c.txt
  gmt grdcut @earth_relief_03s -R$R -Grel_$c.nc
  gmt grdgradient rel_$c.nc -A315 -Nt0.7 -Gint_$c.nc
  gmt grdsample ${c}_ratio.tif -R$R -I6s -nn -Gr6.nc
  gmt grdmath r6.nc 1 SUB ABS 0.03 GT 0 NAN r6.nc MUL = ratio_$c.nc
  gmt grdfilter ${c}_built.tif -R$R -Fg0.6 -D2 -Gbuilt_$c.nc
done
MH=$(gmt mapproject -R$(cat R_kabul.txt) -JM${PW}c -Wh)
GMTV=$(gmt --version); GSV=$(gs --version)

gmt begin $OUT/$NAME pdf,png E300
gmt set FONT_ANNOT_PRIMARY 7p,Helvetica FONT_LABEL 8p,Helvetica MAP_FRAME_TYPE fancy \
        MAP_FRAME_WIDTH 2p MAP_GRID_PEN_PRIMARY 0.25p,gray40 MAP_TICK_LENGTH_PRIMARY 2.5p \
        FORMAT_GEO_MAP ddd:mmF PS_PAGE_COLOR white
for i in 0 1 2 3 4 5; do
  c=${CITIES[$i]}; R=$(cat R_$c.txt); J=M${PW}c
  col=$((i % 3)); row=$((i / 3))
  if [ $i -eq 0 ]; then XY="-Y$(awk -v h=$MH 'BEGIN{print h+6.8}')c"
  elif [ $col -eq 0 ]; then XY="-X-$(awk -v w=$PW -v g=$GAP 'BEGIN{print 2*(w+g)}')c -Y-$(awk -v h=$MH 'BEGIN{print h+0.75}')c"
  else XY="-X$(awk -v w=$PW -v g=$GAP 'BEGIN{print w+g}')c"; fi
  SIDES=WSen
  gmt grdimage int_$c.nc -Cshade.cpt -R$R -J$J $XY
  gmt grdimage ratio_$c.nc -Cratio.cpt -Q -t15
  [ -s ${c}_alluvium.gmt ] && gmt plot ${c}_alluvium.gmt -W0.45p,130/80/30
  [ -s ${c}_fan.gmt ] && gmt plot ${c}_fan.gmt -W0.45p,130/80/30,2_1.5
  gmt grdcontour built_$c.nc -C+1500 -W1.6p,white@25
  gmt grdcontour built_$c.nc -C+1500 -W0.7p,black
  gmt basemap -Bxa0.2f0.1g0.1 -Bya0.2f0.1g0.1 -B$SIDES
  echo ${LON[$c]} ${LAT[$c]} | gmt plot -Sa0.32c -Gwhite -W0.7p,black
  IFS=, read -r _ scen _ pch _ _ _ _ _ rmin rmax rmean _ _ <<< "$(grep "^$c," site_effect.csv)"
  echo "${LET[$i]}) ${TITLE[$c]}" | gmt text -F+cTL+f8.5p,Helvetica-Bold,black+jTL -Dj0.1c/0.1c -Gwhite@10 -C0.06c -N
  awk -v x=${LON[$c]} -v y=${LAT[$c]} -v hx=$HX -v hy=$HY 'BEGIN{printf "> %.4f %.4f 6.5p 4.65c l\n", x-hx+0.012, y-hy+0.075}' > box.txt
  printf "%s; %.0f %% of nodes adjusted\n\nPGA ratio %.2f-%.2f, mean %.2f\n" "$scen" "$pch" "$rmin" "$rmax" "$rmean" >> box.txt
  gmt text box.txt -M -F+f6.5p,Helvetica,black+jBL -Gwhite@10 -C0.06c/0.05c -N
  gmt basemap -Lg$(awk -v x=${LON[$c]} -v y=${LAT[$c]} 'BEGIN{printf "%.3f/%.3f", x+0.2, y+0.2}')+c${LAT[$c]}+w10k+f+lkm \
    -F+gwhite@15+p0.3p,gray40+r0.05c --FONT_ANNOT_PRIMARY=6.5p --FONT_LABEL=6.5p --MAP_SCALE_HEIGHT=4p
done
TW=$(awk -v w=$PW -v g=$GAP 'BEGIN{print 3*w+2*g}')
gmt plot -T -X-$(awk -v w=$PW -v g=$GAP 'BEGIN{print 2*(w+g)}')c
gmt colorbar -Cratio.cpt -Dx$(awk -v w=$TW 'BEGIN{print w/2}')c/-0.9c+w14c/0.26c+h+jTC -L0.06c -S \
  --FONT_ANNOT_PRIMARY=6.5p,Helvetica
gmt text -R0/$TW/0/1 -JX${TW}c/1c -N -F+f8p,Helvetica+jCM <<EOF
$(awk -v w=$TW 'BEGIN{print w/2}') -1.6 PGA ratio: lithology-constrained V@-S30@- / global slope-based V@-S30@- (governing scenario)
EOF
cat > leg.txt <<'EOF'
N 3
S 0.3c - 0.5c - 0.7p,black 0.65c Built-up area, GHSL 2020 (@~\263@~1500 m@+2@+ per 3" cell)
S 0.3c - 0.5c - 0.45p,130/80/30 0.65c Quaternary alluvium (geological map)
S 0.3c - 0.5c - 0.45p,130/80/30,2_1.5 0.65c Quaternary alluvial fan (geological map)
S 0.3c a 0.3c white 0.7p,black 0.65c City centre
EOF
gmt legend leg.txt -Dx$(awk -v w=$TW 'BEGIN{print w/2}')c/-2.05c+w${TW}c+jTC -F+gwhite+p0.4p,gray40 \
  --FONT_ANNOT_PRIMARY=7p,Helvetica -R0/$TW/0/1 -JX${TW}c/1c
gmt text -R0/$TW/0/3 -JX${TW}c/3c -M -N -F+f7p,Helvetica,gray20+jTL -Y-6.0c <<EOF
> 0 2.9 0.3c ${TW}c j
Depicts: change in median PGA when the slope-based USGS Vs30 is constrained by the mapped lithology (alluvium clipped to 180-360 m/s, NEHRP D; alluvial fans to 270-560 m/s; bedrock to @~\263@~360 m/s, NEHRP C or stiffer), for the scenario governing each city (Kabul: Chaman Mw 7.0; Herat and Kandahar: Herat Mw 7.5; Mazar-i-Sharif: Kholm-type Mw 6.5; Kunduz and Jalalabad: Hindu Kush Mw 7.5). Ratios computed at the 30 arc-sec Vs30 nodes with the logic-tree median of openquake.hazardlib 3.26.2 and shown as pixels. Ratios within 0.97-1.03 are not shaded. Box: governing scenario, share of nodes whose Vs30 changes, range of the ratio and its mean over the changed nodes. Built-up outline: GHSL built-up surface @~\263@~1500 m@+2@+ per 3 arc-sec cell after 0.6 km Gaussian smoothing. Mercator projection, windows 0.6@~\260@~ x 0.48@~\260@~. Data: this study; national geological map (USGS/AGS); USGS global Vs30; GHSL GHS-BUILT-S R2023A (European Commission, JRC); SRTM 3 arc-sec (GMT earth_relief) hillshade. Software: GMT $GMTV, GhostScript $GSV, GDAL, Python 3 (openquake.hazardlib, rasterio). Source: authors.
EOF
gmt end
echo done
