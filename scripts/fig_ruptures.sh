#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
ISO=${ISO:-1}
NAME=fig_ruptures
mkdir -p work_rup && cd work_rup

W=60; E=75.2; S=29; N=39.2
R=$W/$E/$S/$N; J=B67.6/34/30.5/37.5/17c
HALO='=~1.2p,white'

curl -s -f -m 180 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
for cls in A B; do
  ogr2ogr -f OGR_GMT qf_$cls.gmt faults.shp -t_srs EPSG:4326 -where "ActvCls LIKE '%($cls)%'"
done
FDSN="https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&minlatitude=$S&maxlatitude=$N&minlongitude=$W&maxlongitude=$E"
curl -s -m 180 "$FDSN&minmagnitude=4&starttime=1973-01-01&endtime=2025-12-31T23:59:59&orderby=time-asc" \
  | awk -F, 'NR>1{print $3,$2,$4,$5}' > seis.txt
NSEIS=$(wc -l < seis.txt)
awk '$3>70' seis.txt > deep.txt

python3 - "$OUT" <<'PY'
import sys, numpy as np, geopandas as gpd, csv
from pyproj import Geod, Proj
geod = Geod(ellps="WGS84")
g = gpd.read_file("faults.shp").to_crs(4326)

def verts(names, lonmin=-999, lonmax=999, cls=None):
    s = g[g.Name.isin(names)]
    if cls: s = s[s.ActvCls.str.contains(cls, regex=False)]
    xy = []
    for geom in s.geometry:
        for part in getattr(geom, "geoms", [geom]):
            xy += [c for c in part.coords if lonmin <= c[0] <= lonmax]
    return np.array(xy)

def fit_line(xy, lon0, lat0):
    """principal axis of the strand vertices in a local azimuthal projection"""
    p = Proj(proj="aeqd", lon_0=lon0, lat_0=lat0)
    x, y = p(xy[:, 0], xy[:, 1])
    X = np.c_[x, y]; c = X.mean(0)
    u = np.linalg.svd(X - c)[2][0]
    t = (X - c) @ u
    return p, c, u, t

def trace_from(p, c, u, t0, t1):
    a = c + t0 * u; b = c + t1 * u
    (lo0, la0), (lo1, la1) = p(*a, inverse=True), p(*b, inverse=True)
    return (lo0, la0), (lo1, la1)

def wc_ss(M):
    return 10 ** (-2.57 + 0.62 * M), 10 ** (-0.76 + 0.27 * M)

def strasser_slab(M):
    return 10 ** (-2.350 + 0.562 * M), 10 ** (-1.058 + 0.356 * M)

rups = []
M = 7.0; L, Wd = wc_ss(M)
xy = verts(["Chaman section E"], cls="(A)")
p, c, u, t = fit_line(xy, 68.5, 33.7)
if (c + t.max() * u)[1] < (c + t.min() * u)[1]: u, t = -u, -t
xyn = xy[t > t.max() - 1.2 * L * 1e3]
p, c, u, t = fit_line(xyn, *xyn.mean(0))
if (c + t.max() * u)[1] < (c + t.min() * u)[1]: u, t = -u, -t
a, b = trace_from(p, c, u, t.max() - L * 1e3, t.max())
rups.append(dict(name="Chaman", M=M, mech="left-lateral", a=a, b=b, dip=90, rake=0, ztop=0.0, W=Wd, L=L,
                 hyp=None))
M = 7.5; L, Wd = wc_ss(M)
xy = verts(["Hari Rud fault system - western section", "Hari Rud fault system - central section"], 61.3, 64.0)
p, c, u, t = fit_line(xy, 62.5, 34.4)
if u[0] < 0: u, t = -u, -t
tc = (np.array(p(62.20, 34.35)) - c) @ u
a, b = trace_from(p, c, u, tc - L * 5e2, tc + L * 5e2)
rups.append(dict(name="Herat", M=M, mech="right-lateral", a=a, b=b, dip=90, rake=180, ztop=0.0, W=Wd, L=L,
                 hyp=None))
M = 7.5; L, Wd = strasser_slab(M)
rups.append(dict(name="Hindu Kush", M=M, mech="reverse, intermediate depth", strike=107.0, dip=71.0, rake=91.0,
                 hyp=(70.3676, 36.5244, 231.0), W=Wd, L=L))
M = 6.5; L, Wd = wc_ss(M)
rups.append(dict(name="Kholm-type", M=M, mech="oblique reverse-strike-slip", strike=241.3, dip=70.2, rake=37.3,
                 hyp=(67.110, 36.709, 21.5), W=Wd, L=L))

def plane_from_centre(lon, lat, z, strike, dip, L, W):
    """rectangle centred on (lon, lat, z); returns top trace (a, b), ztop, zbot, corners"""
    hw = W / 2 * np.cos(np.radians(dip)); dz = W / 2 * np.sin(np.radians(dip))
    ddir = strike + 90
    tlon, tlat, _ = geod.fwd(lon, lat, ddir + 180, hw * 1e3)
    a = geod.fwd(tlon, tlat, strike + 180, L / 2 * 1e3)[:2]
    b = geod.fwd(tlon, tlat, strike, L / 2 * 1e3)[:2]
    return a, b, z - dz, z + dz

rows = []
with open("ruptures.gmt", "w") as fr, open("traces.gmt", "w") as ft, open("meca.txt", "w") as fm:
    for r in rups:
        if r["hyp"] is None:
            (lo0, la0), (lo1, la1) = r["a"], r["b"]
            strike = geod.inv(lo0, la0, lo1, la1)[0] % 360
            r["strike"] = strike; ztop, zbot = 0.0, r["W"]
            clon, clat = (lo0 + lo1) / 2, (la0 + la1) / 2; cz = r["W"] / 2
        else:
            clon, clat, cz = r["hyp"]
            (lo0, la0), (lo1, la1), ztop, zbot = plane_from_centre(clon, clat, cz, r["strike"], r["dip"], r["L"], r["W"])
        horiz = r["W"] * np.cos(np.radians(r["dip"])) * 1e3
        ddir = r["strike"] + 90
        lo2, la2, _ = geod.fwd(lo1, la1, ddir, horiz); lo3, la3, _ = geod.fwd(lo0, la0, ddir, horiz)
        fr.write(f"> -Z{len(rows)}\n{lo0:.4f} {la0:.4f}\n{lo1:.4f} {la1:.4f}\n{lo2:.4f} {la2:.4f}\n{lo3:.4f} {la3:.4f}\n{lo0:.4f} {la0:.4f}\n")
        ft.write(f"> -Z{len(rows)}\n{lo0:.4f} {la0:.4f}\n{lo1:.4f} {la1:.4f}\n")
        fm.write(f"{clon:.4f} {clat:.4f} {cz:.1f} {r['strike']:.0f} {r['dip']:.0f} {r['rake']:.0f} 5 0 0 {r['name']}\n")
        rows.append([r["name"], r["M"], r["mech"], round(r["strike"]), r["dip"], r["rake"], round(ztop, 1), round(zbot, 1),
                     round(r["L"], 1), round(r["W"], 1), round(clon, 3), round(clat, 3), round(cz, 1)])
with open(f"{sys.argv[1]}/ruptures.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["scenario", "Mw", "mechanism", "strike", "dip", "rake", "ztop_km", "zbot_km", "L_km", "W_km",
                "centre_lon", "centre_lat", "centre_depth_km"])
    w.writerows(rows)
for r in rows: print(r)
PY

awk '$1>=68 && $2>=35' deep.txt > deep_hk.txt
gmt blockmedian deep_hk.txt -i0,1,2 -R67.5/75.5/35/39.5 -I0.25 > deep_bm.txt
gmt surface deep_bm.txt -R67.5/75.5/35/39.5 -I0.05 -T0.5 -fg -Gslab_raw.nc
gmt grdfilter slab_raw.nc -Fg70k -D4 -fg -Gslab_sm.nc
awk '{print $1,$2,1}' deep_hk.txt | gmt blockmean -R67.5/75.5/35/39.5 -I0.05 -Ss > cnt.txt
gmt xyz2grd cnt.txt -R67.5/75.5/35/39.5 -I0.05 -fg -Gcnt.nc -di0
gmt grdfilter cnt.nc -Fc60k -D4 -Ni -fg -Gcnt_sm.nc
gmt grdmath cnt_sm.nc 0.044 GE 0 NAN slab_sm.nc MUL = slab.nc

gmt grdcut @earth_relief_15s -R59.5/75.7/28.5/39.7 -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.6 -Gint.nc
[ "$ISO" = 1 ] && gmt grdfilter relief.nc -Fg6k -D4 -I30s -Grelief_c.nc
gmt makecpt -Cgray -T-1/1 -G0.6/0.98 > shade.cpt
cat > depth.cpt <<'CPT'
0   255/215/0   70  255/215/0   ;Crustal, 0-70 km
70  255/140/0   150 255/140/0   ;Upper slab, 70-150 km
150 205/0/0     200 205/0/0     ;Lower slab, 150-200 km
200 125/38/205  700 125/38/205  ;Slab tip, >200 km
B 255/215/0
F 125/38/205
N white
CPT
awk '{m=$4; s=(m<5)?0.06:(m<6)?0.12:(m<7)?0.20:0.30; print $1,$2,$3,s}' seis.txt | sort -k3,3n > seis_cls.txt
gmt coast -EAF -M -Di > af.txt
cat > cities.txt <<'EOF'
69.181 34.519 BL 0.15c/0.06c Kabul
65.710 31.620 TL 0.15c/-0.06c Kandahar
62.200 34.348 TL 0.15c/-0.1c Herat
67.110 36.709 BR -0.15c/0.1c Mazar-i-Sharif
70.450 34.430 TL 0.15c/-0.06c Jalalabad
68.860 36.728 BL 0.15c/0.06c Kunduz
EOF
read MW MH <<< "$(gmt mapproject -R$R -J$J -W)"
GMTV=$(gmt --version); GSV=$(gs --version)
ISOTXT=""; [ "$ISO" = 1 ] && ISOTXT="; topographic contours at 1000 m interval"
COL=("58/95/205" "34/139/34" "0/205/205" "40/40/40")

gmt begin $OUT/$NAME pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p MAP_GRID_PEN_PRIMARY 0.25p,gray45 \
        MAP_TICK_LENGTH_PRIMARY 3p FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel PS_PAGE_COLOR white
gmt grdimage int.nc -Cshade.cpt -R$R -J$J -Y9.6c
gmt coast -EIR,PK,TM,UZ,TJ,CN,IN+gwhite@45
[ "$ISO" = 1 ] && gmt grdcontour relief_c.nc -C1000 -L500/8000 -W0.3p,95/55/25@45 -Q40
gmt coast -Di -W0.25p,gray35 -N1/0.5p,gray30,- -A50 -Bxa2f30mg2 -Bya2f30mg2 -BWSen
gmt plot af.txt -W2.2p,white@20
gmt plot af.txt -W0.9p,gray10
gmt plot qf_B.gmt -W0.3p,gray45
gmt plot qf_A.gmt -W0.7p,gray25
gmt plot seis_cls.txt -Sc -Cdepth.cpt -W0.1p,gray20
gmt grdcontour slab.nc -C50 -A50+f8p,Helvetica,black+gwhite@15 -L90/260 -W0.8p,black,- -Gd6c -Q30
awk 'NR==1{print $1,$2,$3,$4,$5,$6,$7,66.2,33.0}
     NR==2{print $1,$2,$3,$4,$5,$6,$7,62.3,32.55}
     NR==3{print $1,$2,$3,$4,$5,$6,$7,73.3,37.95}
     NR==4{print $1,$2,$3,$4,$5,$6,$7,65.2,37.95}' meca.txt > meca_off.txt
for i in 0 1 2 3; do
  sed -n "$((i+1))p" meca_off.txt | gmt meca -Sa1.05c -G${COL[$i]} -Ewhite -W0.6p,black -A+p0.6p,${COL[$i]}+s0.12c
done
gmt text -F+f8.5p,Helvetica-Bold,black$HALO+j -N <<'EOF'
66.2 32.43 TC Chaman Mw 7.0
62.3 31.98 TC Herat Mw 7.5
73.3 37.38 TC Hindu Kush Mw 7.5
65.2 37.38 TC Kholm-type Mw 6.5
EOF
awk '{print $1,$2}' cities.txt | gmt plot -Ss0.22c -Gwhite -W0.8p,black
for i in 0 1 2 3; do
  gmt convert ruptures.gmt -Q$i | gmt plot -G${COL[$i]}@45 -W0.8p,${COL[$i]}
  gmt convert traces.gmt -Q$i | gmt plot -W2.6p,${COL[$i]}
done
while read lo la j d nm; do
  echo "$lo $la $nm" | gmt text -F+f8.5p,Helvetica-Bold,black$HALO+j$j -D$d
done < cities.txt
gmt basemap -Lg61.9/37.4+c34+w200k+f+l"km" -F+gwhite@15+p0.4p,gray40+r0.1c
gmt basemap -Tdg60.9/38.3+w1.1c+l,,,N
gmt colorbar -Cdepth.cpt -DJBC+w${MW}c/0.3c+o0/0.75c+h -L0.08c
echo "$(awk -v w=$MW 'BEGIN{print w/2}') 0.5 Hypocentral depth" | gmt text -R0/$MW/0/1 -JX${MW}c/1c -Y-2.2c -N -F+f9p,Helvetica+jCM
gmt plot -T -R$R -J$J -Y2.2c
CONT='G 0'; [ "$ISO" = 1 ] && CONT='S 0.3c - 0.55c - 0.3p,95/55/25 0.75c Topographic contour, 1000 m'
cat > leg.txt <<EOF
N 3
S 0.3c r 0.5c/0.32c 58/95/205@45 0.8p,58/95/205 0.75c Rupture plane (surface projection)
S 0.3c - 0.55c - 2.6p,58/95/205 0.75c Rupture top edge
S 0.3c - 0.55c - 0.7p,black,- 0.75c Median depth of seismicity > 70 km
S 0.3c - 0.55c - 0.7p,gray25 0.75c Quaternary fault, class A
S 0.3c - 0.55c - 0.3p,gray45 0.75c Quaternary fault, class B
$CONT
S 0.3c c 0.06c gray55 0.1p,gray20 0.75c Earthquake M 4.0-4.9
S 0.3c c 0.12c gray55 0.1p,gray20 0.75c Earthquake M 5.0-5.9
S 0.3c c 0.20c gray55 0.1p,gray20 0.75c Earthquake M 6.0-6.9
S 0.3c c 0.30c gray55 0.1p,gray20 0.75c Earthquake M @~\263@~ 7.0
S 0.3c s 0.2c white 0.8p,black 0.75c City
EOF
gmt legend leg.txt -DJBC+w${MW}c+o0/2.3c -F+gwhite+p0.4p,gray40 --FONT_ANNOT_PRIMARY=8.5p,Helvetica
gmt text -R0/$MW/0/3 -JX${MW}c/3c -M -N -F+f8p,Helvetica,gray20+jTL -Y-7.2c <<EOF
> 0 2.9 0.35c ${MW}c j
Depicts: surface projections and top edges of the four scenario ruptures with their focal mechanisms (lower-hemisphere, compressional quadrants in the scenario colour): Chaman Mw 7.0 (left-lateral, northern Chaman section E), Herat Mw 7.5 (right-lateral, Hari Rud fault), Hindu Kush Mw 7.5 (repeat of the 26 October 2015 event at 231 km depth) and a Kholm-type Mw 6.5 source beneath Mazar-i-Sharif (mechanism of the 2 November 2025 Kholm earthquake); earthquakes M@~\263@~4, 1973-2025 (n = $NSEIS) in four depth classes and four magnitude classes; contours of the median depth of seismicity deeper than 70 km (50 km interval)${ISOTXT}. Rupture dimensions: Wells and Coppersmith (1994), crustal sources; Strasser et al. (2010), intraslab source. Albers conic equal-area projection. Data: USGS OFR 2007-1103 Quaternary faults; USGS ComCat hypocentres and W-phase moment tensors; SRTM15+ (GMT earth_relief, 15 arc-sec) hillshade; GSHHG; DCW. Software: GMT $GMTV, GhostScript $GSV, Python 3 (geopandas, pyproj, numpy). Source: authors.
EOF
gmt end
echo done
