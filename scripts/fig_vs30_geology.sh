#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
VS30=${VS30:-global_vs30.grd}
GDRIVE_ID=${GDRIVE_ID:-1TLTNVIpkC4p7Mff4VkBo8mUVRfA3Yh6v}
[ -f "$VS30" ] && VS30=$(realpath "$VS30")
mkdir -p work_vs && cd work_vs

W=60; E=75.2; S=29; N=39.2
R=$W/$E/$S/$N; J=B67.6/34/30.5/37.5/8.2c
HALO='=~1.1p,white'

if [ -f "$VS30" ]; then
  gmt grdcut "$VS30" -R59.5/75.7/28.5/39.7 -Gvs30.nc
else
  U=https://apps.usgs.gov/shakemap_geodata/vs30/global_vs30.grd
  GDAL_SKIP=netCDF gdal_translate -q -srcwin 28740 10140 1945 1345 /vsicurl/$U raw.tif
  python3 - <<'PY'
import rasterio, numpy as np
from rasterio.transform import from_origin
a = rasterio.open('raw.tif').read(1)[::-1]
d = 1/120
with rasterio.open('vs30.tif', 'w', driver='GTiff', height=a.shape[0], width=a.shape[1], count=1,
                   dtype='float32', crs='EPSG:4326',
                   transform=from_origin(59.5 - d/2, 39.7 + d/2, d, d), nodata=np.nan) as o:
    o.write(a.astype('float32'), 1)
PY
  gdal_translate -q -of netCDF vs30.tif vs30.nc
fi
gmt grdinfo vs30.nc -C

curl -s -f -m 600 -L -o geol.zip "https://drive.usercontent.google.com/download?id=$GDRIVE_ID&export=download&confirm=t"
unzip -q -o geol.zip '*/Geologic age/*'
GEO=$(ls */"Geologic age"/naf_disk/shapes/geology/glg*.shp | head -1)
python3 - "$GEO" <<'PY'
import sys, re, geopandas as gpd
g = gpd.read_file(sys.argv[1])
def cls(s):
    s = (s or '').strip()
    if s.lower().startswith('w'): return 'water'
    if re.match(r'^Q\d*(ac|t)$', s): return 'fan'
    if re.match(r'^Q\d*(a|lap|lac|sm|e|loe)$', s): return 'alluvium'
    return 'bedrock'
g['CLS'] = g['UNIT_SYMBO'].map(cls)
q = g[g['UNIT_SYMBO'].str.startswith('Q', na=False)]
print(q.groupby(['UNIT_SYMBO', 'CLS']).size().to_string())
print(g['CLS'].value_counts().to_string())
for c in ['alluvium', 'fan', 'bedrock', 'water']:
    g[g.CLS == c][['CLS', 'geometry']].to_file(f'geo_{c}.gpkg', driver='GPKG')
PY
for c in alluvium fan bedrock water; do ogr2ogr -f OGR_GMT geo_$c.gmt geo_$c.gpkg; done

python3 - <<'PY'
import numpy as np, rasterio, geopandas as gpd
from rasterio.features import rasterize
src = rasterio.open('vs30.nc')
v = src.read(1).astype(float)
for k, c in enumerate(['alluvium', 'fan', 'bedrock'], 1):
    g = gpd.read_file(f'geo_{c}.gpkg')
    m = rasterize(((geom, 1) for geom in g.geometry), out_shape=v.shape, transform=src.transform) == 1
    x = v[m & np.isfinite(v)]
    print(f'{c:9s} n={x.size:7d} median={np.median(x):6.0f} IQR={np.percentile(x,25):5.0f}-{np.percentile(x,75):5.0f}'
          f'  <360: {100*np.mean(x<360):4.1f}%  360-760: {100*np.mean((x>=360)&(x<760)):4.1f}%  >=760: {100*np.mean(x>=760):4.1f}%')
PY

gmt grdcut @earth_relief_15s -R59.5/75.7/28.5/39.7 -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.6 -Gint.nc
REG=$(gmt grdinfo vs30.nc -C | awk '{print ($12==1) ? "-r" : ""}')
gmt grdsample int.nc -Rvs30.nc -I30s $REG -Gint_vs.nc
gmt grdmath int_vs.nc 0.45 MUL = int_soft.nc
gmt makecpt -Croma -T150/900/10 -D > vs30.cpt
gmt makecpt -Cgray -T-1/1 -G0.55/0.98 > shade.cpt
gmt coast -EAF -M -Di > af.txt

cat > cities.txt <<'EOF'
69.181 34.519 MR -0.13c/0.06c Kabul
65.710 31.620 ML 0.13c/0 Kandahar
62.200 34.348 ML 0.13c/0 Herat
67.110 36.709 MR -0.13c/0.06c Mazar-i-Sharif
70.450 34.430 ML 0.13c/-0.08c Jalalabad
68.860 36.728 ML 0.13c/0.06c Kunduz
EOF
read MW MH <<< "$(gmt mapproject -R$R -J$J -W)"
GAP=0.6; TW=$(awk -v w=$MW -v g=$GAP 'BEGIN{print 2*w+g}'); DX=$(awk -v w=$MW -v g=$GAP 'BEGIN{print w+g}')
GMTV=$(gmt --version); GSV=$(gs --version); GDALV=$(gdalinfo --version | awk '{print $2}' | tr -d ,)

cities() {
  awk '{print $1,$2}' cities.txt | gmt plot -Ss0.17c -Gred3 -W0.5p,white
  while read lo la j d nm; do
    echo "$lo $la $nm" | gmt text -F+f8p,Helvetica-Bold,black$HALO+j$j -D$d
  done < cities.txt
}
frame() {
  gmt coast -Di -W0.25p,gray30 -N1/0.4p,gray25,- -A50 -Bxa4f1g2 -Bya2f1g2 -B$1
  gmt plot af.txt -W1.6p,white@20
  gmt plot af.txt -W0.7p,gray10
}

gmt begin $OUT/fig_vs30_geology pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 2p MAP_GRID_PEN_PRIMARY 0.25p,gray45 \
        MAP_TICK_LENGTH_PRIMARY 2.5p FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel PS_PAGE_COLOR white

gmt grdimage vs30.nc -Iint_soft.nc -Cvs30.cpt -R$R -J$J -Y7c
gmt coast -EIR,PK,TM,UZ,TJ,CN,IN+gwhite@42
frame WSen
cities
echo "60.45 39.0 (a)" | gmt text -F+f11p,Helvetica-Bold,black$HALO+jTL
gmt basemap -Tdg74.1/30.2+w0.8c+l,,,N
gmt colorbar -Cvs30.cpt -DJBC+w${MW}c/0.25c+o0/0.75c+h+e -Bxa100f20+l"@%2%V@%%@-S30@- (m/s)"
gmt plot -R150/900/0/1 -JX${MW}c/0.25c -Y-1.0c -W0.9p,black <<'EOF'
>
180 -0.2
180 1.2
>
360 -0.2
360 1.2
>
760 -0.2
760 1.2
EOF
gmt text -R150/900/0/1 -JX${MW}c/0.25c -F+f8p,Helvetica-Bold,black$HALO+jCM -N <<'EOF'
270 0.5 D
560 0.5 C
830 0.5 B
EOF
gmt plot -T -Y1.0c

gmt grdimage int.nc -Cshade.cpt -R$R -J$J -X${DX}c
gmt plot geo_bedrock.gmt  -G165/150/135 -t30
gmt plot geo_fan.gmt      -G235/150/60  -t15
gmt plot geo_alluvium.gmt -G250/225/120 -t15
gmt plot geo_water.gmt    -G110/170/220
gmt coast -EIR,PK,TM,UZ,TJ,CN,IN+gwhite@42
frame wSEn
cities
echo "60.45 39.0 (b)" | gmt text -F+f11p,Helvetica-Bold,black$HALO+jTL
gmt basemap -Lg73.0/30.0+c34+w200k -F+gwhite@15+p0.4p,gray40+r0.1c --FONT_ANNOT_PRIMARY=8p --FONT_LABEL=8p
cat > leg.txt <<'EOF'
N 1
S 0.3c r 0.5c/0.3c 250/225/120 0.4p,gray40 0.75c Alluvium (alluvial, playa, eolian sand, loess)
S 0.3c r 0.5c/0.3c 235/150/60 0.4p,gray40 0.75c Alluvial fans, colluvium and till
S 0.3c r 0.5c/0.3c 165/150/135 0.4p,gray40 0.75c Bedrock (pre-Quaternary and volcanic rocks)
S 0.3c r 0.5c/0.3c 110/170/220 0.4p,gray40 0.75c Water body
S 0.3c s 0.17c red3 0.5p,white 0.75c Scenario city
EOF
gmt legend leg.txt -DJBC+w${MW}c+o0/0.75c -F+gwhite+p0.4p,gray40 --FONT_ANNOT_PRIMARY=8p,Helvetica

gmt text -R0/$TW/0/3 -JX${TW}c/3c -M -N -F+f8p,Helvetica,gray20+jTL -X-${DX}c -Y-5.9c <<EOF
> 0 2.9 0.35c ${TW}c j
Depicts: (a) time-averaged shear-wave velocity of the upper 30 m (V@-S30@-) from the topographic-slope proxy, with NEHRP site-class boundaries at 180, 360 and 760 m/s (classes D, C, B); (b) geological map units grouped into Quaternary alluvium, alluvial fans with colluvium and till, and bedrock; the six scenario cities. Albers conic equal-area projection (central meridian 67.6@.E, standard parallels 30.5@.N and 37.5@.N). Data: USGS global V@-S30@- (30 arc-sec); USGS national geological map (2006); SRTM15+ (GMT earth_relief, 15 arc-sec) hillshade; GSHHG; DCW. Software: GMT $GMTV, GDAL $GDALV, GhostScript $GSV, Python 3 (geopandas, rasterio). Source: authors.
EOF
gmt end
echo done
