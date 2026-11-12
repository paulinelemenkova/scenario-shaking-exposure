#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
mkdir -p work_sa_c && cd work_sa_c

W=60; E=75.2; S=29; N=39.2
R=$W/$E/$S/$N; J=B67.6/34/30.5/37.5/17c
HALO='=~1.2p,white'

curl -s -f -m 180 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
for cls in A B C; do
  ogr2ogr -f OGR_GMT qf_R_$cls.gmt  faults.shp -t_srs EPSG:4326 -where "Type='R' AND ActvCls LIKE '%($cls)%'"
  ogr2ogr -f OGR_GMT qf_SS_$cls.gmt faults.shp -t_srs EPSG:4326 -where "Type IN ('RL','LL','SS') AND ActvCls LIKE '%($cls)%'"
  ogr2ogr -f OGR_GMT qf_N_$cls.gmt  faults.shp -t_srs EPSG:4326 -where "Type='N' AND ActvCls LIKE '%($cls)%'"
done

curl -s -m 120 -L -o pb_b.txt http://peterbird.name/oldFTP/PB2002/PB2002_boundaries.dig.txt
awk -F, '/^\*\*\*/{next} /^[^ +-]/{print ">"; next} {gsub(/ /,""); printf "%.4f %.4f\n",$1,$2}' pb_b.txt > pb_bound.gmt

FDSN="https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&minlatitude=$S&maxlatitude=$N&minlongitude=$W&maxlongitude=$E"
curl -s -m 120 "$FDSN&minmagnitude=4&mindepth=70&starttime=1900-01-01&endtime=2025-12-31T23:59:59" \
  | awk -F, 'NR>1{print $3,$2,$4,$5}' > deep.txt
NDEEP=$(wc -l < deep.txt)
: > ev.txt
while read yr t0 t1 mmin; do
  curl -s -m 60 "$FDSN&minmagnitude=$mmin&starttime=$t0&endtime=$t1&orderby=magnitude" \
    | awk -F, -v y=$yr 'NR>1{printf "%s %s %s %s %s\n",$3,$2,$4,$5,y}' >> ev.txt
done <<'EOF'
2015 2015-10-26 2015-10-27 7.0
2022 2022-06-21 2022-06-23 5.8
2023 2023-10-07 2023-10-16 6.2
2025 2025-08-31 2025-09-01 5.8
2025 2025-11-02 2025-11-04 6.0
EOF
cat ev.txt

gmt grdcut @earth_relief_15s -R59.5/75.7/28.5/39.7 -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.6 -Gint.nc
gmt makecpt -Cetopo1 -G0/8500 -T0/7000/100 > topo.cpt
gmt coast -EAF -M -Di > af.txt

cat > cities.txt <<'EOF'
69.181 34.519 BL Kabul
65.710 31.620 TL Kandahar
62.200 34.348 BL Herat
67.110 36.709 BC Mazar-i-Sharif
70.450 34.430 TL Jalalabad
68.860 36.728 BL Kunduz
EOF
read MW MH <<< "$(gmt mapproject -R$R -J$J -W)"
GMTV=$(gmt --version); GSV=$(gs --version)

gmt begin $OUT/fig_studyarea pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica FONT_TITLE 9p,Helvetica \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p MAP_GRID_PEN_PRIMARY 0.25p,gray45 \
        MAP_TICK_LENGTH_PRIMARY 3p FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel PS_PAGE_COLOR white
gmt grdimage relief.nc -Iint.nc -Ctopo.cpt -R$R -J$J -Y9.6c
gmt coast -EIR,PK,TM,UZ,TJ,CN,IN+gwhite@42
gmt coast -Di -W0.25p,gray30 -N1/0.5p,gray25,- -I1/0.6p,30/100/190 -I2/0.35p,30/100/190 -A50 \
          -Bxa2f30mg2 -Bya2f30mg2 -BWSen
gmt plot af.txt -W2.4p,white@20
gmt plot af.txt -W0.9p,gray10

for t in "R 150/20/20" "SS 20/20/20" "N 0/110/60"; do
  set -- $t
  gmt plot qf_$1_C.gmt -W0.2p,$2
  gmt plot qf_$1_B.gmt -W0.45p,$2
  gmt plot qf_$1_A.gmt -W1.0p,$2
done

gmt plot deep.txt -Sc0.08c -G110/40/160@20 -W0.1p,white
gmt plot ev.txt -i0,1 -Sa0.45c -Gyellow -W0.6p,black

awk '{print $1,$2}' cities.txt | gmt plot -Ss0.24c -Gred3 -W0.7p,white
while read lo la j nm; do
  case $j in BL) d=0.18c/0.06c;; TL) d=0.18c/-0.06c;; BC) d=0/0.2c;; esac
  echo "$lo $la $nm" | gmt text -F+f9p,Helvetica-Bold,black$HALO+j$j -D$d
done < cities.txt

awk 'NR==1{printf "%s %s BL 0.25c/0.12c %s M %.1f, %d km deep\n",$1,$2,$5,$4,$3}
     NR==2{printf "%s %s TL 0.25c/-0.1c %s M %.1f\n",$1,$2,$5,$4}
     NR==3{printf "%s %s BC 0/0.3c %s M %.1f sequence (four events)\n",$1,$2,$5,$4}
     NR==7{printf "%s %s BL 0.25c/0.1c %s M %.1f\n",$1,$2,$5,$4}
     NR==8{printf "%s %s TC 0/-0.3c %s M %.1f\n",$1,$2,$5,$4}' ev.txt > evlab.txt
while read lo la j d txt; do
  echo "$lo $la $txt" | gmt text -F+f8p,Helvetica,black$HALO+j$j -D$d
done < evlab.txt

gmt text -F+f+a+j -N <<EOF
67.80 32.75 8p,Helvetica-Oblique,20/20/20$HALO 60 CM Chaman fault
64.10 34.62 8p,Helvetica-Oblique,20/20/20$HALO 0 CM Herat fault
69.20 35.80 8p,Helvetica-Oblique,150/20/20$HALO 28 CM Hindu Kush
66.80 34.95 8p,Helvetica-Oblique,black$HALO 0 CM Koh-i-Baba
63.40 35.55 8p,Helvetica-Oblique,black$HALO 0 CM Paropamisus
73.60 38.45 8p,Helvetica-Oblique,black$HALO 0 CM Pamir
65.10 30.55 8p,Helvetica-Oblique,92/51/23$HALO 0 CM Registan
62.95 30.20 8p,Helvetica-Oblique,92/51/23$HALO 0 CM Dasht-e Margo
61.10 31.15 8p,Helvetica-Oblique,black$HALO 0 CM Sistan Basin
66.00 37.55 8p,Helvetica-Oblique,30/90/180$HALO -18 CM Amu Darya
63.90 31.45 8p,Helvetica-Oblique,30/90/180$HALO -20 CM Helmand R.
68.45 34.62 8p,Helvetica-Oblique,30/90/180$HALO 0 CM Kabul R.
EOF
gmt text -F+f9p,Helvetica,gray15$HALO+j -N <<'EOF'
60.90 35.70 CM Iran
68.80 30.40 CM Pakistan
63.10 38.35 CM Turkmenistan
66.70 38.70 CM Uzbekistan
69.90 38.55 CM Tajikistan
74.80 37.75 CM China
EOF

gmt basemap -Lg62.4/37.15+c34+w200k+f+l"km" -F+gwhite@15+p0.4p,gray40+r0.1c
gmt basemap -Tdg60.9/38.25+w1.1c+l,,,N

gmt inset begin -DjBR+w4.2c/3.1c+o0.35c/0.35c -F+gwhite+p0.5p,gray30
  gmt coast -R38/102/8/50 -JM4.2c -Gwheat2 -Slightsteelblue1 -Dc -A5000 -EAF+gred3 -W0.15p,gray40 -B0
  gmt plot pb_bound.gmt -W0.7p,purple3
  printf "$W $S\n$E $S\n$E $N\n$W $N\n$W $S\n" | gmt plot -W0.8p,black -A
  gmt text -F+f7p,Helvetica-Bold,gray15+j -N <<'EOF'
78 20 CM India
52 33 CM Iran
70 45 CM Eurasia
EOF
gmt inset end

gmt colorbar -Ctopo.cpt -DJBC+w${MW}c/0.3c+o0/0.75c+h+ef -Bxa1000f100+l"Elevation (m)" -I0.4
cat > leg.txt <<EOF
N 3
S 0.3c - 0.55c - 1.0p,150/20/20 0.75c Thrust / reverse fault
S 0.3c - 0.55c - 1.0p,20/20/20 0.75c Strike-slip fault
S 0.3c - 0.55c - 1.0p,0/110/60 0.75c Normal fault
S 0.3c - 0.55c - 1.0p,gray30 0.75c Activity A: > 10 mm/yr
S 0.3c - 0.55c - 0.45p,gray30 0.75c Activity B: 1-10 mm/yr
S 0.3c - 0.55c - 0.2p,gray30 0.75c Activity C: indeterminate
S 0.3c c 0.12c 110/40/160 0.1p,white 0.75c Earthquake, depth > 70 km
S 0.3c a 0.35c yellow 0.6p,black 0.75c Reference earthquake
S 0.3c s 0.22c red3 0.6p,white 0.75c Scenario city
EOF
gmt legend leg.txt -DJBC+w${MW}c+o0/2.3c -F+gwhite+p0.4p,gray40 --FONT_ANNOT_PRIMARY=8.5p,Helvetica
gmt text -R0/$MW/0/2 -JX${MW}c/2c -M -N -F+f8p,Helvetica,gray20+jTL -Y-5.6c <<EOF
> 0 1.45 0.35c ${MW}c j
Depicts: relief of the study domain and surroundings; Quaternary faults by kinematics and activity class; earthquakes deeper than 70 km (M@~\263@~4, 1900-2025, n = $NDEEP) outlining the Hindu Kush intermediate-depth nest; mainshocks (USGS ComCat) of the 2015 Hindu Kush, 2022 Paktika, 2023 Herat and 2025 Kunar and Kholm earthquakes; the six scenario cities. Albers conic equal-area projection (central meridian 67.6@.E, standard parallels 30.5@.N and 37.5@.N); inset with PB2002 plate boundaries. Data: SRTM15+ (GMT earth_relief, 15 arc-sec) with the ETOPO1 colour palette; USGS OFR 2007-1103 (Ruleman et al. 2007); USGS ComCat; Bird (2003) PB2002; GSHHG; DCW. Software: GMT $GMTV, GDAL, GhostScript $GSV. Source: authors.
EOF
gmt end
echo done
