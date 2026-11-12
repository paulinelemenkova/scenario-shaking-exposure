#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
NAME=fig_validation_maps
mkdir -p work_val && cd work_val
EVS=(us7000hj3u us6000ldpg us7000qsvj us6000rl31)
fetch() { if [[ $1 == http* ]]; then curl -s -f -o $2 $1/$2 || true; else cp $1/$2 . 2>/dev/null || true; fi; }
for e in "${EVS[@]}"; do for s in sm.tif model.tif res.tif rupture.gmt stations.txt; do fetch $VAL ${e}_$s; done; done
fetch $VAL validation_stats.csv
TITLE=("2022 Paktika M@-w@- 6.0" "2023 Herat M@-w@- 6.3" "2025 Kunar M@-w@- 6.0" "2025 Kholm M@-w@- 6.2")
ROWLAB=("USGS ShakeMap" "This study (model)" "ln(ShakeMap / model)")
LET=(a b c d e f g h i j k l)
PW=3.55; GAP=1.05; VG=0.45

gmt makecpt -Cturbo -T-0.30103/2/0.05 -Z > pga.cpt
cat > pga_ann.txt <<'EOF'
-0.30103 a 0.5
0 a 1
0.30103 a 2
0.69897 a 5
1 a 10
1.30103 a 20
1.69897 a 50
2 a 100
EOF
gmt makecpt -Cvik -T-1.2/1.2/0.2 -D > res.cpt
cat > cities.txt <<'EOF'
69.181 34.519 Kabul
62.200 34.348 Herat
67.110 36.709 Mazar
70.450 34.430 Jalalabad
68.860 36.728 Kunduz
69.920 33.340 Khost
71.140 34.870 Asadabad
65.740 36.660 Shibirghan
70.130 33.130 Urgun
EOF
GMTV=$(gmt --version); GSV=$(gs --version)
for e in "${EVS[@]}"; do
  IFS=, read -r W E S N <<< "$(awk -F, -v e=$e 'NR>1 && $1==e{print $(NF-5)","$(NF-4)","$(NF-3)","$(NF-2)}' validation_stats.csv)"
  R=$W/$E/$S/$N; echo $R > R_$e.txt
  gmt grdcut @earth_relief_30s -R$R -Grel_$e.nc
  gmt grdgradient rel_$e.nc -A315 -Nt0.5 -Gi_$e.nc
  for k in sm model; do
    gmt grdmath ${e}_$k.tif 100 MUL LOG10 = l_${e}_$k.nc
    gmt grdsample l_${e}_$k.nc -R$R -I30s -rg -Gs_${e}_$k.nc
  done
  gmt grdsample ${e}_res.tif -R$R -I30s -rg -Gs_${e}_res.nc
  gmt grdsample i_$e.nc -R$R -I30s -rg -Gis_$e.nc
done
MH=$(for e in "${EVS[@]}"; do gmt mapproject -R$(cat R_$e.txt) -JM${PW}c -Wh; done | sort -g | tail -1)

gmt begin $OUT/$NAME pdf,png E300
gmt set FONT_ANNOT_PRIMARY 6.5p,Helvetica FONT_LABEL 7.5p,Helvetica MAP_FRAME_TYPE fancy \
        MAP_FRAME_WIDTH 2p MAP_GRID_PEN_PRIMARY 0.25p,gray40 MAP_TICK_LENGTH_PRIMARY 2.5p \
        FORMAT_GEO_MAP ddd:mmF PS_PAGE_COLOR white MAP_ANNOT_OFFSET_PRIMARY 2p
n=0
for row in 0 1 2; do
  for col in 0 1 2 3; do
    e=${EVS[$col]}; R=$(cat R_$e.txt); J=M${PW}c
    if [ $n -eq 0 ]; then XY="-X1.0c -Y$(awk -v h=$MH 'BEGIN{print 2*h+7.2}')c"
    elif [ $col -eq 0 ]; then XY="-X-$(awk -v w=$PW -v g=$GAP 'BEGIN{print 3*(w+g)}')c -Y-$(awk -v h=$MH -v v=$VG 'BEGIN{print h+v}')c"
    else XY="-X$(awk -v w=$PW -v g=$GAP 'BEGIN{print w+g}')c"; fi
    S=$( [ $row -eq 2 ] && echo WSen || echo Wsen )
    gmt basemap -R$R -J$J -B+n $XY
    case $row in
      0) gmt grdimage s_${e}_sm.nc -Cpga.cpt -Iis_$e.nc ;;
      1) gmt grdimage s_${e}_model.nc -Cpga.cpt -Iis_$e.nc ;;
      2) gmt grdimage s_${e}_res.nc -Cres.cpt -Iis_$e.nc ;;
    esac
    gmt coast -Di -W0.25p,gray25 -N1/0.6p,black,- -A100 -Bxa1f0.5g0.5 -Bya1f0.5g0.5 -B$S
    gmt plot ${e}_rupture.gmt -Gwhite@40 -W0.7p,black
    IFS=, read -r elon elat <<< "$(awk -F, -v e=$e 'NR>1 && $1==e{print $(NF-1)","$NF}' validation_stats.csv)"
    echo $elon $elat | gmt plot -Sa0.3c -Gwhite -W0.6p,black
    if [ $row -eq 0 ]; then
      awk '$3=="macroseismic"{print $1,$2}' ${e}_stations.txt | gmt plot -Sc0.08c -Gwhite -W0.3p,black
      awk '$3=="seismic"{print $1,$2}' ${e}_stations.txt | gmt plot -St0.22c -Gmagenta -W0.5p,black
    fi
    IFS=/ read -r w0 e0 s0 n0 <<< "$R"
    awk -v w=$w0 -v e=$e0 -v s=$s0 -v n=$n0 '$1>w+0.05 && $1<e-0.45 && $2>s+0.05 && $2<n-0.12' cities.txt > cty.txt
    gmt plot cty.txt -Ss0.12c -Gwhite -W0.5p,black -i0,1
    gmt text cty.txt -F+f5.8p,Helvetica,black=~0.8p,white+jBL -D0.06c/0.04c
    echo "${LET[$n]})" | gmt text -F+cTL+f7.5p,Helvetica-Bold,black+jTL -Dj0.06c/0.06c -Gwhite@10 -C0.04c -N
    if [ $row -eq 0 ]; then
      echo "${TITLE[$col]}" | gmt text -F+cTC+f8p,Helvetica-Bold,black+jBC -D0/0.12c -N
    fi
    if [ $row -eq 2 ]; then
      IFS=, read -r nn b sd <<< "$(awk -F, -v e=$e 'NR>1 && $1==e{print $(NF-12)","$(NF-11)","$(NF-10)}' validation_stats.csv)"
      printf "bias %+.2f, @~s@~ %.2f\n" "$b" "$sd" | gmt text -F+cBR+f6.3p,Helvetica,black+jBR -Dj0.06c/0.06c -Gwhite@10 -C0.04c -N
    fi
    if [ $col -eq 0 ]; then
      echo "${ROWLAB[$row]}" | gmt text -F+cLM+f7.5p,Helvetica-Bold,black+jBC+a90 -D-0.75c/0 -N
    fi
    n=$((n+1))
  done
done
TW=$(awk -v w=$PW -v g=$GAP 'BEGIN{print 4*w+3*g}')
gmt plot -T -X-$(awk -v w=$PW -v g=$GAP 'BEGIN{print 3*(w+g)}')c
gmt colorbar -Cpga.cpt -Dx$(awk -v w=$TW 'BEGIN{print w*0.25}')c/-0.95c+w7.2c/0.25c+h+jTC+e -Bxcpga_ann.txt \
  -Bx+l"PGA (% g), rows 1-2"
gmt colorbar -Cres.cpt -Dx$(awk -v w=$TW 'BEGIN{print w*0.75}')c/-0.95c+w7.2c/0.25c+h+jTC+e -Bxa0.4f0.2 \
  -Bx+l"ln(PGA@-ShakeMap@- / PGA@-model@-), row 3"
cat > leg.txt <<'EOF'
N 5
S 0.25c r 0.4c/0.25c white@40 0.7p,black 0.5c Rupture plane
S 0.25c a 0.25c white 0.6p,black 0.5c Epicentre
S 0.25c t 0.2c magenta 0.5p,black 0.5c Seismic station
S 0.25c c 0.08c white 0.3p,black 0.5c Macroseismic (DYFI)
S 0.25c s 0.12c white 0.5p,black 0.5c City or town
EOF
gmt legend leg.txt -Dx$(awk -v w=$TW 'BEGIN{print w/2}')c/-2.25c+w${TW}c+jTC -F+gwhite+p0.4p,gray40 \
  --FONT_ANNOT_PRIMARY=6.8p,Helvetica -R0/$TW/0/1 -JX${TW}c/1c
gmt text -R0/$TW/0/3 -JX${TW}c/3c -M -N -F+f6.8p,Helvetica,gray20+jTL -Y-6.2c <<EOF
> 0 2.9 0.28c ${TW}c j
Depicts: USGS ShakeMap (version 4) PGA and the PGA predicted by the method of this study for four recent earthquakes, in 3@~\260@~ x 3@~\260@~ windows centred on the epicentre. Model: logic-tree median of Boore et al. (2014), Cauzzi et al. (2015) and Akkar et al. (2014) (weights 0.4/0.3/0.3), computed with openquake.hazardlib 3.26.2 at every third ShakeMap node (~0.05@~\260@~) with the ShakeMap's Vs30; ruptures: ShakeMap finite faults for Paktika and Kunar, W-phase nodal planes with Wells and Coppersmith (1994) dimensions for Herat and Kholm. Residual bias and standard deviation (ln units) over nodes with R@-JB@- @~\243@~ 200 km and PGA@-ShakeMap@- @~\263@~ 0.5 % g. Seismic stations: instrumental records used by ShakeMap; macroseismic: Did You Feel It? intensities. Mercator projection. Data: USGS ShakeMap, ComCat and W-phase moment tensors; SRTM30+ (GMT earth_relief, 30 arc-sec) hillshade; GSHHG. Software: GMT $GMTV, GhostScript $GSV, Python 3 (openquake.hazardlib, rasterio). Source: authors.
EOF
gmt end
echo done
