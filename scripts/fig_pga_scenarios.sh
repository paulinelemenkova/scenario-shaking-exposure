#!/usr/bin/env bash
set -e
IN=${IN:-data}
OUT=${OUT:-$(pwd)}
NAME=fig_pga_scenarios
mkdir -p work_pga && cd work_pga
if [ -n "$BASE" ]; then
  mkdir -p in; IN=$(pwd)/in
  for f in chaman_pga.tif herat_pga.tif hindukush_pga.tif kholm_pga.tif rup.gmt top.gmt; do curl -s -o in/$f $BASE/$f; done
fi

W=60; E=75.2; S=29; N=39.2
R=$W/$E/$S/$N; PW=8.45; J=B67.6/34/30.5/37.5/${PW}c
HALO='=~1p,white'

gmt grdcut @earth_relief_30s -R59.5/75.7/28.5/39.7 -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.6 -Gint.nc
gmt makecpt -Cgray -T-1/1 -G0.62/0.98 > shade.cpt
LO=(2.76 6.2 11.5 21.5 40.1); HI=(6.2 11.5 21.5 40.1 74.7)
LAB=("IV, light" "V, moderate" "VI, strong" "VII, very strong" "VIII, severe")
CC=($(gmt makecpt -Cturbo -T0/5/1 -G0.22/0.95 -N | awk '!/^#/{print $2}'))
for k in 0 1 2 3 4; do echo "${LO[$k]} ${CC[$k]} ${HI[$k]} ${CC[$k]} ;${LAB[$k]}"; done > pga.cpt
cat pga.cpt
gmt coast -EAF -M -Di > af.txt
cat > cities.txt <<'EOF'
69.181 34.519 Kabul
65.710 31.620 Kandahar
62.200 34.348 Herat
67.110 36.709 Mazar-i-Sharif
70.450 34.430 Jalalabad
68.860 36.728 Kunduz
EOF
read MW MH <<< "$(gmt mapproject -R$R -J$J -W)"
GMTV=$(gmt --version); GSV=$(gs --version)

printf "%s c\n" 6.2 11.5 21.5 40.1 > cont.txt
TAGS=(chaman herat hindukush kholm)
TITLE=("Chaman M@-w@- 7.0" "Herat M@-w@- 7.5" "Hindu Kush M@-w@- 7.5, 231 km" "Kholm-type M@-w@- 6.5")
LET=(a b c d)
for i in 0 1 2 3; do
  t=${TAGS[$i]}
  gmt grdmath $IN/${t}_pga.tif 100 MUL LOG10 = lg.nc
  gmt grdsample lg.nc -R59.5/75.7/28.5/39.7 -I30s -Glg30.nc
  gmt grdmath 10 lg30.nc POW DUP 2.76 GE 0 NAN MUL = ${t}.nc
  [ -f int30.nc ] || gmt grdsample int.nc -R${t}.nc -Gint30.nc
  gmt grd2xyz $IN/${t}_pga.tif -s -o2 | gmt info -C -o1 | awk '{printf "%.2f\n",$1}' > max_${t}.txt
done

gmt begin $OUT/$NAME pdf,png E300
gmt set FONT_ANNOT_PRIMARY 7.5p,Helvetica FONT_LABEL 8.5p,Helvetica \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 2.5p MAP_GRID_PEN_PRIMARY 0.25p,gray45 \
        MAP_TICK_LENGTH_PRIMARY 2.5p FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel PS_PAGE_COLOR white
for i in 0 1 2 3; do
  t=${TAGS[$i]}
  col=$((i % 2)); row=$((i / 2))
  if [ $i -eq 0 ]; then gmt basemap -R$R -J$J -B+n -Y$(awk -v h=$MH 'BEGIN{print h+7.0}')c
  elif [ $i -eq 1 ]; then gmt basemap -R$R -J$J -B+n -X$(awk -v w=$MW 'BEGIN{print w+0.45}')c
  elif [ $i -eq 2 ]; then gmt basemap -R$R -J$J -B+n -X-$(awk -v w=$MW 'BEGIN{print w+0.45}')c -Y-$(awk -v h=$MH 'BEGIN{print h+0.45}')c
  else gmt basemap -R$R -J$J -B+n -X$(awk -v w=$MW 'BEGIN{print w+0.45}')c; fi
  FR="-Bxa4f1g2 -Bya2f1g2"
  SIDES=$( [ $col -eq 0 ] && echo W || echo w )$( [ $row -eq 1 ] && echo S || echo s )en
  gmt grdimage int.nc -Cshade.cpt
  gmt grdimage ${t}.nc -Cpga.cpt -Q -t12
  gmt grdcontour ${t}.nc -Ccont.txt -W0.25p,gray15@35 -Q40
  gmt coast -EIR,PK,TM,UZ,TJ,CN,IN+gwhite@50
  gmt coast -Di -W0.25p,gray35 -N1/0.4p,gray30,- -A50 $FR -B$SIDES
  gmt plot af.txt -W1.6p,white@20
  gmt plot af.txt -W0.7p,gray10
  gmt convert $IN/rup.gmt -Q$i | gmt plot -Gwhite@60 -W0.8p,black
  gmt convert $IN/top.gmt -Q$i | gmt plot -W1.0p,black
  awk '{print $1,$2}' cities.txt | gmt plot -Ss0.16c -Gwhite -W0.6p,black
  awk '$3=="Kabul"||$3=="Herat"||$3=="Kandahar"{print $1,$2,$3}
       $3=="Mazar-i-Sharif"{print $1,$2,"Mazar"}' cities.txt \
    | gmt text -F+f7p,Helvetica-Bold,black$HALO+jBL -D0.1c/0.05c
  echo "${LET[$i]}) ${TITLE[$i]}" | gmt text -F+cTL+f8.5p,Helvetica-Bold,black+jTL -Dj0.12c/0.12c -Gwhite@15 -C0.06c -N
  echo "PGA@-max@- $(cat max_${t}.txt) g" | gmt text -F+cBR+f7.5p,Helvetica,black+jBR -Dj0.12c/0.12c -Gwhite@15 -C0.06c -N
done
TW=$(awk -v w=$MW 'BEGIN{print 2*w+0.45}')
gmt plot -T -X-$(awk -v w=$MW 'BEGIN{print w+0.45}')c
gmt colorbar -Cpga.cpt -Dx$(awk -v w=$TW 'BEGIN{print w/2}')c/-0.95c+w14c/0.28c+h+jTC -L0.06c -S \
  --FONT_ANNOT_PRIMARY=7p,Helvetica
gmt text -R0/$TW/0/1 -JX${TW}c/1c -N -F+f8.5p,Helvetica+jCM <<EOF
$(awk -v w=$TW 'BEGIN{print w/2}') -1.75 Median PGA as instrumental intensity, MMI (class limits 2.8, 6.2, 11.5, 21.5, 40.1 and 74.7 % g; Worden et al., 2012)
EOF
gmt text -R0/$TW/0/3 -JX${TW}c/3c -M -N -F+f7.5p,Helvetica,gray20+jTL -Y-5.0c <<EOF
> 0 2.9 0.32c ${TW}c j
Depicts: median PGA of the four scenario earthquakes: logic-tree weighted geometric mean of ground-motion models (crustal: Boore et al., 2014, 0.4; Cauzzi et al., 2015, 0.3; Akkar et al., 2014, 0.3; intraslab: Parker et al., 2022, 0.5; Kuehn et al., 2023, 0.5) computed with openquake.hazardlib 3.26.2 on a 0.05@~\260@~ grid with site Vs30 from the USGS global slope-based grid; resampled to 30 arc-sec in log space; PGA below 2.8 % g (MMI II-III, weak) is not shaded; thin lines: class limits. White polygon: surface projection of the rupture plane; thick line: top edge. Albers conic equal-area projection. Data: USGS global Vs30; SRTM30+ (GMT earth_relief, 30 arc-sec) hillshade; GSHHG; DCW. Software: GMT $GMTV, GhostScript $GSV, Python 3 (openquake.hazardlib, rasterio, numpy). Source: authors.
EOF
gmt end
echo done
