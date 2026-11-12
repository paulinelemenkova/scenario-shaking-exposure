#!/usr/bin/env bash
set -e
IN=${IN:-data}
OUT=${OUT:-$(pwd)}
NAME=fig_exposure_map
mkdir -p work_exp && cd work_exp
if [ -n "$BASE" ]; then
  mkdir -p in; IN=$(pwd)/in
  for f in chaman_pga.tif herat_pga.tif hindukush_pga.tif kholm_pga.tif rup.gmt top.gmt; do curl -s -o in/$f $BASE/$f; done
fi

G=https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E2025_GLOBE_R2023A_4326_30ss/V1-0/tiles
P=GHS_POP_E2025_GLOBE_R2023A_4326_30ss_V1_0
for t in R5_C25 R5_C26 R6_C25 R6_C26 R7_C25 R7_C26; do
  [ -f ${P}_$t.tif ] || { curl -s -o $t.zip $G/${P}_$t.zip; unzip -o -q $t.zip ${P}_$t.tif; }
done
gdalbuildvrt -q pop.vrt ${P}_R*.tif
gdal_translate -q -of netCDF -projwin 59.5 39.7 75.7 28.5 pop.vrt popraw.nc
gmt grdmath popraw.nc 0 MAX = pop.nc
RP=$(gmt grdinfo pop.nc -Ir); IP=$(gmt grdinfo pop.nc -I)
gmt grdmath pop.nc Y COSD 0.858640 MUL DIV DUP 10 GE 0 NAN MUL = dens.nc
gmt coast -EAF -M -Di > af.txt
gmt grdmask af.txt $RP $IP -r -NNaN/1/1 -Gafmask.nc
gmt grdmath pop.nc afmask.nc MUL = popaf.nc

TAGS=(chaman herat hindukush kholm)
python3 - "$IN" "${TAGS[@]}" <<'PY'
"""Bilinear resampling of log10 PGA onto the GHSL cells and people per ShakeMap MMI class."""
import csv
import sys

import numpy as np
import rasterio
from rasterio.warp import Resampling, reproject

src_dir, tags = sys.argv[1], sys.argv[2:]
lo = [0.0, 2.76, 6.2, 11.5, 21.5, 40.1, 74.7]
lab = ["<IV", "IV", "V", "VI", "VII", "VIII", "IX+"]
with rasterio.open("popaf.nc") as p:
    pop = np.nan_to_num(p.read(1).astype(float))
    tr, shape = p.transform, p.shape
tot = pop.sum()
rows = []
for t in tags:
    with rasterio.open(f"{src_dir}/{t}_pga.tif") as s:
        lg = np.log10(s.read(1).astype(np.float64) * 100.0)
        dst = np.full(shape, np.nan)
        reproject(lg, dst, src_transform=s.transform, src_crs="EPSG:4326", dst_transform=tr,
                  dst_crs="EPSG:4326", resampling=Resampling.bilinear, src_nodata=np.nan, dst_nodata=np.nan)
    with rasterio.open(f"lg30_{t}.tif", "w", driver="GTiff", height=shape[0], width=shape[1], count=1,
                       dtype="float32", crs="EPSG:4326", transform=tr, nodata=np.nan) as o:
        o.write(dst.astype("float32"), 1)
    pga = np.nan_to_num(10**dst, nan=0.0)
    cls = np.digitize(pga, lo[1:])
    e = [pop[cls == k].sum() for k in range(len(lab))]
    rows.append([t, *[f"{v:.0f}" for v in e]])
    ge5, ge6 = sum(e[2:]), sum(e[3:])

    def fmt(v: float) -> str:
        if v < 500:
            return "0"
        return f"{v / 1e6:.2f} M" if v >= 1e6 else f"{v / 1e3:.0f} k"

    with open(f"lab_{t}.txt", "w") as f:
        f.write(f"{fmt(ge6)}|{fmt(ge5)}\n")
    print(t, "PGAmax", round(float(np.nanmax(10**dst)), 1), "MMI>=V", round(ge5), "MMI>=VI", round(ge6),
          [round(v) for v in e])
with open("exposure_national.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["scenario", *lab])
    w.writerows(rows)
    w.writerow(["total_domain", f"{tot:.0f}", *[""] * (len(lab) - 1)])
print("total", round(tot))
PY
cp exposure_national.csv $OUT/

W=60; E=75.2; S=29; N=39.2
R=$W/$E/$S/$N; PW=15.0; J=B67.6/34/30.5/37.5/${PW}c
HALO='=~1.2p,white'
gmt grdcut @earth_relief_30s -R59.5/75.7/28.5/39.7 -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.6 -Gint.nc
gmt makecpt -Cgray -T-1/1 -G0.80/1.0 > shade.cpt
DL=(10 50 250 1000 5000); DH=(50 250 1000 5000 100000)
DLAB=("10-50" "50-250" "250-1000" "1000-5000" "@~\263@~5000")
GR=(195 150 105 60 15)
for k in 0 1 2 3 4; do g=${GR[$k]}; c="$g/$((g-3))/$((g-8))"; echo "${DL[$k]} $c ${DH[$k]} $c ;${DLAB[$k]}"; done > dens.cpt
cat > cities.txt <<'EOF'
69.181 34.519 Kabul BL 69.62 34.88
65.710 31.620 Kandahar BL - -
62.200 34.348 Herat TC 62.200 33.42
67.110 36.709 Mazar-i-Sharif BC 66.85 37.62
70.450 34.430 Jalalabad TL - -
68.860 36.728 Kunduz TC 68.55 35.88
EOF
COL=("213/94/0" "0/114/178" "0/158/115" "204/121/167")
NAMES=("Chaman M@-w@- 7.0" "Herat M@-w@- 7.5" "Hindu Kush M@-w@- 7.5" "Kholm-type M@-w@- 6.5")
read MW MH <<< "$(gmt mapproject -R$R -J$J -W)"
GMTV=$(gmt --version); GSV=$(gs --version)
TOTM=$(awk -F, '$1=="total_domain"{printf "%.1f", $2/1e6}' exposure_national.csv)
echo "0.7924 C" > c5.txt
echo "1.0607 C" > c6.txt

gmt begin $OUT/$NAME pdf,png E300
gmt set FONT_ANNOT_PRIMARY 7.5p,Helvetica FONT_LABEL 8.5p,Helvetica \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 2.5p MAP_GRID_PEN_PRIMARY 0.25p,gray45 \
        MAP_TICK_LENGTH_PRIMARY 2.5p FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel PS_PAGE_COLOR white
gmt basemap -R$R -J$J -B+n -Y9.0c
gmt grdimage int.nc -Cshade.cpt
gmt grdimage dens.nc -Cdens.cpt -Q
gmt coast -EIR,PK,TM,UZ,TJ,CN,IN,KG+gwhite@45
gmt coast -Di -W0.25p,gray35 -N1/0.4p,gray30,- -A50 -Bxa2f1g2 -Bya2f1g2 -BWSen
gmt plot af.txt -W1.6p,white@20
gmt plot af.txt -W0.7p,gray10
for i in 0 1 2 3; do
  t=${TAGS[$i]}; c=${COL[$i]}
  gmt grdcontour lg30_$t.tif -Cc5.txt -W0.8p,$c,4_2:0 -Q60
  gmt grdcontour lg30_$t.tif -Cc6.txt -W1.0p,$c -Q60
done
gmt plot $IN/rup.gmt -Gwhite@50 -W0.6p,black
gmt plot $IN/top.gmt -W1.0p,black
awk '$5!="-"{printf ">\n%s %s\n%s %s\n", $1, $2, $5, $6}' cities.txt | gmt plot -W0.5p,gray15
awk '{print $1,$2}' cities.txt | gmt plot -Ss0.18c -Gwhite -W0.6p,black
while read -r x y n j lx ly; do
  case $j in BL) d="0.12c/0.06c";; TL) d="0.12c/-0.06c";; BR) d="-0.12c/0.06c";; *) d="0c/0c";; esac
  if [ "$lx" = "-" ]; then echo "$x $y $n" > lab.txt; else echo "$lx $ly $n" > lab.txt; d="0c/0c"; fi
  gmt text lab.txt -F+f7.5p,Helvetica-Bold,black$HALO+j$j -D$d
done < cities.txt
gmt basemap -LjBR+c34+w200k+f+o0.35c/0.45c+l"km" -F+gwhite@20+p0.3p,gray40 --FONT_ANNOT_PRIMARY=7p,Helvetica --FONT_LABEL=7p,Helvetica
gmt basemap -TdjTL+w0.9c+f2+l,,,N+o0.6c/0.45c --FONT_TITLE=8p,Helvetica-Bold

gmt plot -T -Y-0.95c
cat > leg.txt <<EOF
N 2
EOF
for i in 0 1 2 3; do
  IFS="|" read -r E6 E5 < lab_${TAGS[$i]}.txt
  echo "S 0.45c - 0.75c - 1.0p,${COL[$i]} 1.0c ${NAMES[$i]}: @~\263@~VI ${E6}; @~\263@~V ${E5}" >> leg.txt
done
gmt legend leg.txt -Dx0c/0c+w${PW}c+jTL+l1.35 -F+gwhite+p0.3p,gray40 -C0.12c/0.08c --FONT_ANNOT_PRIMARY=7.5p,Helvetica
gmt plot -T -Y-0.75c
cat > leg2.txt <<EOF
N 2
S 0.45c - 0.75c - 1.0p,gray20 1.0c MMI VI isoline (PGA 11.5 % g, onset of damage)
S 0.45c - 0.75c - 0.8p,gray20,4_2:0 1.0c MMI V isoline (PGA 6.2 % g)
EOF
gmt legend leg2.txt -Dx0c/0c+w${PW}c+jTL+l1.35 -F+gwhite+p0.3p,gray40 -C0.12c/0.08c --FONT_ANNOT_PRIMARY=7.5p,Helvetica
gmt colorbar -Cdens.cpt -Dx$(awk -v w=$PW 'BEGIN{print w/2}')c/-1.0c+w10c/0.26c+h+jTC -L0.06c -S \
  --FONT_ANNOT_PRIMARY=7p,Helvetica
gmt text -R0/$PW/0/1 -JX${PW}c/1c -N -F+f8p,Helvetica+jCM <<EOF
$(awk -v w=$PW 'BEGIN{print w/2}') -1.85 Population density, people km@+-2@+ (GHSL 2025, 30 arc-sec; < 10 km@+-2@+ not shaded)
EOF
gmt text -R0/$PW/0/3 -JX${PW}c/3c -M -N -F+f7.5p,Helvetica,gray20+jTL -Y-5.1c <<EOF
> 0 2.9 0.32c ${PW}c j
Depicts: population exposed to the four scenario earthquakes. Grey: GHSL population density (GHS-POP R2023A, epoch 2025, 30 arc-sec), neighbouring countries faded. Coloured lines: isolines of the logic-tree median PGA of each scenario (openquake.hazardlib 3.26.2, USGS global Vs30, resampled to 30 arc-sec in log space) at the ShakeMap limits of MMI VI (11.5 % g, solid) and MMI V (6.2 % g, dashed) (Worden et al., 2012); the Hindu Kush event (231 km depth) does not reach MMI VI. Legend: people in the domain (total ${TOTM} million in GHSL 2025) inside each isoline. White polygon: surface projection of the rupture plane; thick line: top edge. Albers conic equal-area projection. Data: GHSL (European Commission, JRC); USGS global Vs30; SRTM30+ (GMT earth_relief, 30 arc-sec) hillshade; GSHHG; DCW. Software: GMT $GMTV, GhostScript $GSV, GDAL, Python 3 (openquake.hazardlib, rasterio, numpy). Source: authors.
EOF
gmt end
echo done
