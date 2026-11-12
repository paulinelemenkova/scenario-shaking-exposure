import csv
import os
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from pyproj import Geod
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score
from sklearn.model_selection import GroupKFold

IN = "data/val"
OUT = "data/res"
os.makedirs(OUT, exist_ok=True)
geod = Geod(ellps="WGS84")
rng = np.random.default_rng(7)
EVS = ["us7000hj3u", "us6000ldpg", "us7000qsvj", "us6000rl31"]

df = pd.read_csv(f"{IN}/validation_nodes.csv")
parts = []
for ev in EVS:
    d = df[df.event == ev].copy()

    with rasterio.open(f"{IN}/{ev}_cls.tif") as src:
        r, c = rasterio.transform.rowcol(src.transform, d.lon.values, d.lat.values)
        a = src.read(1)
        ok = (
            (np.array(r) >= 0)
            & (np.array(r) < a.shape[0])
            & (np.array(c) >= 0)
            & (np.array(c) < a.shape[1])
        )
        cls = np.full(len(d), 255)
        cls[ok] = a[np.array(r)[ok], np.array(c)[ok]]
    d["lith"] = cls

    xy = np.loadtxt(f"{IN}/{ev}_rupture.gmt", comments=">")
    strike = geod.inv(xy[0, 0], xy[0, 1], xy[1, 0], xy[1, 1])[0]
    clon, clat = xy[:4, 0].mean(), xy[:4, 1].mean()
    az = geod.inv(
        np.full(len(d), clon), np.full(len(d), clat), d.lon.values, d.lat.values
    )[0]
    d["az_strike"] = np.abs(np.sin(np.radians(az - strike)))

    st: np.ndarray = np.empty((0, 2))
    if os.path.getsize(f"{IN}/{ev}_stations.txt") > 0:
        st = np.loadtxt(f"{IN}/{ev}_stations.txt", usecols=(0, 1), ndmin=2)
    if len(st):
        dmin = np.full(len(d), np.inf)
        for x, y in st:
            _, _, dist = geod.inv(
                d.lon.values, d.lat.values, np.full(len(d), x), np.full(len(d), y)
            )
            dmin = np.minimum(dmin, dist / 1e3)
        d["d_obs"] = np.minimum(dmin, 999.0)
    else:
        d["d_obs"] = 999.0
    d["res_within"] = d.res_ln - d.res_ln.mean()
    parts.append(d)
df = pd.concat(parts, ignore_index=True)

rows = []
rb = np.geomspace(1, 200, 13)
vb = np.arange(150, 951, 75)
for ev in EVS:
    d = df[df.event == ev]
    for key, bins in (("rjb_km", rb), ("vs30", vb)):
        idx = np.digitize(d[key], bins) - 1
        for k in range(len(bins) - 1):
            s = d.res_ln[idx == k]
            if len(s) >= 15:
                centre = (
                    np.sqrt(bins[k] * bins[k + 1])
                    if key == "rjb_km"
                    else (bins[k] + bins[k + 1]) / 2
                )
                rows.append(
                    [
                        ev,
                        key,
                        f"{centre:.3f}",
                        len(s),
                        f"{s.mean():.4f}",
                        f"{s.std(ddof=1):.4f}",
                    ]
                )
with open(f"{OUT}/binned.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["event", "x", "centre", "n", "mean", "sd"])
    w.writerows(rows)

pts = pd.concat(
    [
        df[df.event == ev].sample(
            n=min(120, int((df.event == ev).sum())), random_state=7
        )
        for ev in EVS
    ]
)
pts[["event", "rjb_km", "vs30", "res_ln"]].to_csv(
    f"{OUT}/points.csv", index=False, float_format="%.2f"
)

m = df.lith.isin([1, 2, 3])
X = pd.DataFrame(
    {
        "log10 R_JB": np.log10(df.rjb_km[m].clip(lower=1)),
        "Vs30": df.vs30[m],
        "Lithology": df.lith[m],
        "Azimuth to strike": df.az_strike[m],
        "Distance to observation": np.log10(df.d_obs[m].clip(lower=1)),
    }
)
y = df.res_within[m].values
blocks = (
    np.floor(df.lon[m] / 0.5).astype(int) * 1000 + np.floor(df.lat[m] / 0.5).astype(int)
).values
pred = np.zeros_like(y)
imp: dict[str, list[float]] = {str(k): [] for k in X.columns}
for tr, te in GroupKFold(n_splits=5).split(X, y, blocks):
    rf = RandomForestRegressor(
        n_estimators=300, min_samples_leaf=20, n_jobs=-1, random_state=7
    )
    rf.fit(X.iloc[tr], y[tr])
    pred[te] = rf.predict(X.iloc[te])
    pi = permutation_importance(
        rf, X.iloc[te], y[te], n_repeats=10, random_state=7, n_jobs=-1
    )
    for col, val in zip(X.columns, pi.importances_mean, strict=True):
        imp[str(col)].append(float(val))
r2 = r2_score(y, pred)
with open(f"{OUT}/importance.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["feature", "mean", "sd"])
    for name, vals in sorted(imp.items(), key=lambda kv: -np.mean(kv[1])):
        w.writerow([name, f"{np.mean(vals):.4f}", f"{np.std(vals, ddof=1):.4f}"])
with open(f"{OUT}/rf_score.txt", "w") as f:
    f.write(f"n={len(y)} r2_cv={r2:.3f}\n")
print("n", len(y), "R2 (spatial-block CV)", round(r2, 3))
print(Path(f"{OUT}/importance.csv").read_text())
print(df.groupby(["event", "lith"]).res_ln.agg(["count", "mean"]).round(3))
