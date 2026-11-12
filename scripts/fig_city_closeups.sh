#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
NAME=fig_city_closeups
mkdir -p work_city && cd work_city
CITIES=(kabul herat mazar kunduz jalalabad kandahar)
fetch() { if [[ $1 == http* ]]; then curl -s -f -o $2 $1/$2 || true; else cp $1/$2 . 2>/dev/null || true; fi; }
for c in "${CITIES[@]}"; do
  fetch $ENV ${c}_env.tif; fetch $ENV ${c}_gov.tif
  fetch $PREP ${c}_built.tif; fetch $PREP ${c}_alluvium.gmt; fetch $PREP ${c}_fan.gmt
done
fetch $ENV city_pga.csv; fetch $RUP rup.gmt; fetch $RUP top.gmt

declare -A LON LAT TITLE
LON=([kabul]=69.181 [herat]=62.200 [mazar]=67.110 [kunduz]=68.860 [jalalabad]=70.450 [kandahar]=65.710)
LAT=([kabul]=34.519 [herat]=34.348 [mazar]=36.709 [kunduz]=36.728 [jalalabad]=34.430 [kandahar]=31.620)
TITLE=([kabul]=Kabul [herat]=Herat [mazar]=Mazar-i-Sharif [kunduz]=Kunduz [jalalabad]=Jalalabad [kandahar]=Kandahar)
LET=(a b c d e f)
PW=5.05; GAP=1.35; HX=0.30; HY=0.24

LO=(0.3 2.76 6.2 11.5 21.5 40.1); HI=(2.76 6.2 11.5 21.5 40.1 74.7)
LAB=("II-III, weak" "IV, light" "V, moderate" "VI, strong" "VII, very strong" "VIII, severe")
CC=("214/226/236" $(gmt makecpt -Cturbo -T0/5/1 -G0.22/0.95 -N | awk '!/^#/{print $2}'))
for k in 0 1 2 3 4 5; do echo "${LO[$k]} ${CC[$k]} ${HI[$k]} ${CC[$k]} ;${LAB[$k]}"; done > pga.cpt
printf "%s c\n" 2.76 6.2 11.5 21.5 40.1 > cont.txt
printf "%s c\n" 0.5 1.5 2.5 > govc.txt
gmt makecpt -Cgray -T-1/1 -G0.55/1 > shade.cpt
mmi() { awk -v p=$1 'BEGIN{split("0.3 2.76 6.2 11.5 21.5 40.1 74.7",b," ");
  split("I II-III IV V VI VII VIII",m," "); r="IX+"; for(i=1;i<=7;i++) if(p<b[i]){r=m[i]; break}; print r}'; }

for c in "${CITIES[@]}"; do
  lo=${LON[$c]}; la=${LAT[$c]}
  R=$(awk -v x=$lo -v y=$la -v hx=$HX -v hy=$HY 'BEGIN{printf "%.4f/%.4f/%.4f/%.4f", x-hx, x+hx, y-hy, y+hy}')
  echo $R > R_$c.txt
  gmt grdcut @earth_relief_03s -R$R -Grel_$c.nc
  gmt grdgradient rel_$c.nc -A315 -Nt0.7 -Gint_$c.nc
  gmt grdmath ${c}_env.tif 100 MUL LOG10 = lg.nc
  gmt grdsample lg.nc -R$R -I6s -Glg6.nc
  gmt grdmath 10 lg6.nc POW DUP 0.3 GE 0 NAN MUL = pga_$c.nc
  gmt grdsample ${c}_gov.tif -R$R -I6s -nn -Ggov_$c.nc
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
  gmt grdimage pga_$c.nc -Cpga.cpt -Q -t30
  gmt grdcontour pga_$c.nc -Ccont.txt -W0.3p,gray20@30
  [ -s ${c}_alluvium.gmt ] && gmt plot ${c}_alluvium.gmt -W0.45p,130/80/30
  [ -s ${c}_fan.gmt ] && gmt plot ${c}_fan.gmt -W0.45p,130/80/30,2_1.5
  gmt grdcontour built_$c.nc -C+1500 -W1.6p,white@25
  gmt grdcontour built_$c.nc -C+1500 -W0.7p,black
  gmt plot rup.gmt -Gwhite@65 -W0.7p,black
  gmt plot top.gmt -W1.2p,black
  gmt basemap -Bxa0.2f0.1g0.1 -Bya0.2f0.1g0.1 -B$SIDES
  echo ${LON[$c]} ${LAT[$c]} | gmt plot -Sa0.32c -Gwhite -W0.7p,black
  IFS=, read -r _ scen mw pga _ _ _ vs <<< "$(grep "^$c," city_pga.csv)"
  P=$(awk -v g=$pga 'BEGIN{printf "%.1f", 100*g}')
  echo "${LET[$i]}) ${TITLE[$c]}" | gmt text -F+cTL+f8.5p,Helvetica-Bold,black+jTL -Dj0.1c/0.1c -Gwhite@10 -C0.06c -N
  awk -v x=${LON[$c]} -v y=${LAT[$c]} -v hx=$HX -v hy=$HY 'BEGIN{printf "> %.4f %.4f 6.5p 4.65c l\n", x-hx+0.012, y-hy+0.075}' > box.txt
  printf "Governing: %s M@-w@- %s\n\nCentre: PGA %s %% g, MMI %s; V@-S30@- %s m/s\n" "$scen" "$mw" "$P" "$(mmi $P)" "$vs" >> box.txt
  gmt text box.txt -M -F+f6.5p,Helvetica,black+jBL -Gwhite@10 -C0.06c/0.05c -N
  gmt basemap -Lg$(awk -v x=${LON[$c]} -v y=${LAT[$c]} 'BEGIN{printf "%.3f/%.3f", x+0.2, y+0.2}')+c${LAT[$c]}+w10k+f+lkm \
    -F+gwhite@15+p0.3p,gray40+r0.05c --FONT_ANNOT_PRIMARY=6.5p --FONT_LABEL=6.5p --MAP_SCALE_HEIGHT=4p
done
TW=$(awk -v w=$PW -v g=$GAP 'BEGIN{print 3*w+2*g}')
gmt plot -T -X-$(awk -v w=$PW -v g=$GAP 'BEGIN{print 2*(w+g)}')c
gmt colorbar -Cpga.cpt -Dx$(awk -v w=$TW 'BEGIN{print w/2}')c/-0.9c+w15c/0.26c+h+jTC -L0.06c -S \
  --FONT_ANNOT_PRIMARY=6.5p,Helvetica
gmt text -R0/$TW/0/1 -JX${TW}c/1c -N -F+f8p,Helvetica+jCM <<EOF
$(awk -v w=$TW 'BEGIN{print w/2}') -1.6 Envelope of the median PGA over the four scenarios as instrumental intensity, MMI (Worden et al., 2012)
EOF
cat > leg.txt <<'EOF'
N 3
S 0.3c - 0.5c - 0.7p,black 0.65c Built-up area, GHSL 2020 (@~\263@~1500 m@+2@+ per 3" cell)
S 0.3c - 0.5c - 0.45p,130/80/30 0.65c Quaternary alluvium (geological map)
S 0.3c - 0.5c - 0.45p,130/80/30,2_1.5 0.65c Quaternary alluvial fan (geological map)
S 0.3c r 0.5c/0.3c white@65 0.7p,black 0.65c Rupture plane (thick line: top edge)
S 0.3c a 0.3c white 0.7p,black 0.65c City centre
EOF
gmt legend leg.txt -Dx$(awk -v w=$TW 'BEGIN{print w/2}')c/-2.05c+w${TW}c+jTC -F+gwhite+p0.4p,gray40 \
  --FONT_ANNOT_PRIMARY=7p,Helvetica -R0/$TW/0/1 -JX${TW}c/1c
gmt text -R0/$TW/0/3 -JX${TW}c/3c -M -N -F+f7p,Helvetica,gray20+jTL -Y-6.0c <<EOF
> 0 2.9 0.3c ${TW}c j
Depicts: city-scale shaking for the six largest cities of the domain. Colour: maximum over the four scenarios (Chaman Mw 7.0, Herat Mw 7.5, Hindu Kush Mw 7.5, Kholm-type Mw 6.5) of the logic-tree median PGA, computed at the 30 arc-sec nodes of the USGS global Vs30 grid with openquake.hazardlib 3.26.2 and resampled to 6 arc-sec in log space; thin grey lines: class limits. Box: scenario governing shaking at the city centre, with its PGA, MMI class and Vs30. Built-up outline: GHSL built-up surface @~\263@~1500 m@+2@+ per 3 arc-sec cell after 0.6 km Gaussian smoothing. Mercator projection, windows 0.6@~\260@~ x 0.48@~\260@~. Data: this study; GHSL GHS-BUILT-S R2023A, epoch 2020, 3 arc-sec (European Commission, JRC); national geological map (USGS/AGS); USGS global Vs30; SRTM 3 arc-sec (GMT earth_relief) hillshade. Software: GMT $GMTV, GhostScript $GSV, GDAL, Python 3 (openquake.hazardlib, geopandas, rasterio). Source: authors.
EOF
gmt end
echo done
