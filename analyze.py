#!/usr/bin/env python3
"""
BiteCast cloud analyzer — Florida Keys
--------------------------------------
Runs on a schedule (GitHub Actions), pulls free satellite + bathymetry grids,
detects temperature/colour FRONTS and bottom STRUCTURE, fuses them into a
fish-suitability map, layers a per-species SPAWNING-WINDOW model, then applies
an FWC / national-marine-sanctuary CONSERVATION filter before it ever suggests
a spot. It pushes a summary to your phone (ntfy) and writes hotspots.geojson
(which the web app can load).

Data (all free, no key):
  * SST + Chlorophyll come from a LADDER of NOAA ERDDAP datasets (see SST_SOURCES /
                 CHL_SOURCES). Each run tries them in order of quality, reads the date of
                 the newest image, REFUSES anything older than its max_age_days, and logs
                 which one it used. (The MODIS-Aqua feeds this tool was born on stopped
                 updating in 2022 — a fixed dataset id is a trap; a ladder with a freshness
                 check is not.) Daily 4 km chlorophyll images are stacked into a 7-day
                 median composite to fill cloud holes.
  * Bathymetry   NOAA ERDDAP  srtm15plus       (~500 m — folds in NOAA's Coastal Relief
                 Model near the coast; ~3.7x finer than the old ETOPO)
  * Currents     NOAA CoastWatch  noaacwBLENDEDNRTcurrentsDaily  (altimetry geostrophic,
                 near-real-time, 0.25°) — gives the Florida-Current edge, eddies & convergence
  * Sargassum    USF/AOML  noaa_aoml_atlantic_oceanwatch_AFAI_7D  (AFAI, 7-day, ~0.015°)
                 — floating-algae index for weed lines (mahi/tripletail/bait)
  * Structure    FWC artificial reefs (4,400+ pts) + NOAA AWOIS wrecks/obstructions
                 (ArcGIS FeatureServer, GeoJSON) — exact reef/wreck coordinates.
                 DISPLAY/PLANNING ONLY — neither agency verifies every point. NOT for navigation.

NOTHING here overrides the law. Regulations change constantly — this tool is a
planning aid, not legal advice. Always confirm current rules in the FWC
"Fish Rules" app and at myfwc.com before fishing, and respect every closure.
"""

import os, re, sys, json, math, datetime, warnings, urllib.parse
import numpy as np
import requests

# ERDDAP servers (different datasets live on different nodes)
ERDDAP_PFEG = "https://coastwatch.pfeg.noaa.gov/erddap"   # SST, chl, bathymetry
ERDDAP_CW   = "https://coastwatch.noaa.gov/erddap"        # blended near-real-time currents
ERDDAP_AOML = "https://cwcgom.aoml.noaa.gov/erddap"       # USF/AOML sargassum (AFAI)

# ArcGIS REST layers for structure (queried with an envelope, returned as GeoJSON)
FWC_REEFS  = "https://gis.myfwc.com/mapping/rest/services/Open_Data/Artificial_Reef_Locations_in_Florida/MapServer/12"
NOAA_WRECKS = "https://maps.nccs.nasa.gov/mapping/rest/services/hifld_open/transportation_water/FeatureServer/39"

# ============================== CONFIG ==============================
CONFIG = {
    # Whole Florida Keys box (extend lonmin to -83.2 to include the Dry Tortugas)
    "box": {"latmin": 24.30, "latmax": 25.60, "lonmin": -82.40, "lonmax": -80.00},
    "stride": 1,                       # ERDDAP index stride (raise to 2 to thin the grid)
    "ntfy_topic": os.environ.get("NTFY_TOPIC", "bitecast-CHANGE-ME-7h3xq9"),
    "target_species": "mutton",        # see SPECIES below
    "home": {"name": "Islamorada", "lat": 24.92, "lon": -80.63},  # for "near you" ranking
    "top_n": 6,
    # fusion weights (fish-finding map) — auto-renormalised over whatever layers load
    "w_front": 0.34,        # coincident temp + colour break
    "w_structure": 0.14,    # bottom slope / ledges (from bathymetry)
    "w_reef": 0.12,         # proximity to known wrecks / artificial reefs
    "w_depth": 0.10,        # species depth band
    "w_current": 0.14,      # current edge / shear (the Florida-Current "wall")
    "w_converge": 0.10,     # surface convergence (bait + weed accumulate here)
    "w_weed": 0.06,         # sargassum / floating-algae signal
    "min_edge_score": 0.20,
    # extra signals (set False to skip; the score renormalises automatically)
    "use_currents": True,
    "use_sargassum": True,
    "use_structure": True,             # FWC reefs + NOAA wrecks
    # data sources for the extra signals
    "current_dataset": "noaacwBLENDEDNRTcurrentsDaily",  # 0.25° NRT geostrophic, u_current/v_current (m/s)
    "sargassum_dataset": "noaa_aoml_atlantic_oceanwatch_AFAI_7D",  # USF/AOML 7-day AFAI, var "AFAI"
    "bathy_dataset": "srtm15plus",     # ~500 m; var "z". Swap to a 90 m CRM volume for even finer detail.
    "bathy_var": "z",
    "bathy_stride": 2,                 # srtm15plus is dense; stride 2 ≈ ~900 m (still ~2x finer than ETOPO)
    "structure_radius_nm": 0.75,       # how close to a reef/wreck counts as "on structure"
}

# ======================= SST / CHLOROPHYLL LADDERS =======================
# Tried top to bottom; the first dataset that answers AND whose newest image is no older
# than max_age_days wins. "days" > 1 pulls that many daily slices and takes the per-pixel
# median (fills cloud holes). Verified on the NOAA catalog 2026-09-26 unless marked.
SST_SOURCES = [
    {"label": "MUR 1 km (gap-free analysis)",          "server": ERDDAP_PFEG, "dataset": "jplMURSST41",
     "var": "analysed_sst", "days": 1, "max_age_days": 4},     # was 5 weeks behind on 2026-09-26 — kept on top for when it catches up
    {"label": "Geo-Polar blended 5 km (gap-free)",     "server": ERDDAP_PFEG, "dataset": "nesdisGeoPolarSSTN5NRT",
     "var": "analysed_sst", "days": 1, "max_age_days": 4},
    {"label": "Coral Reef Watch 5 km",                 "server": ERDDAP_PFEG, "dataset": "NOAA_DHW",
     "var": "CRW_SST", "days": 1, "max_age_days": 7},
]
CHL_SOURCES = [
    {"label": "VIIRS NOAA-20 4 km, 7-day composite",   "server": ERDDAP_CW,   "dataset": "noaacwN20VIIRSchlaDaily",
     "var": "chlor_a", "days": 7, "max_age_days": 5},          # origin server — id from the catalog, unverified from the sandbox
    {"label": "VIIRS S-NPP 4 km, 7-day composite",     "server": ERDDAP_CW,   "dataset": "noaacwNPPVIIRSchlaDaily",
     "var": "chlor_a", "days": 7, "max_age_days": 5},          # origin server — unverified from the sandbox
    {"label": "VIIRS NOAA-20 4 km, 7-day composite (mirror)", "server": ERDDAP_PFEG, "dataset": "nesdisVHNnoaa20chlaDaily",
     "var": "chlor_a", "days": 7, "max_age_days": 5},
    {"label": "VIIRS S-NPP 4 km, 7-day composite (mirror)",   "server": ERDDAP_PFEG, "dataset": "nesdisVHNchlaDaily",
     "var": "chlor_a", "days": 7, "max_age_days": 5},
    {"label": "VIIRS gap-filled 9 km (DINEOF)",        "server": ERDDAP_PFEG, "dataset": "nesdisVHNnoaaSNPPnoaa20NRTchlaGapfilledDaily",
     "var": "chlor_a", "days": 1, "max_age_days": 5},          # coarse but reliably current — the safety net
]

# ============================ SPECIES =============================
# General Keys guidance (spawning season / lunar tendency / temp window °F / depth band ft).
# Generalised from the literature — NOT precise, and many of these aggregate on
# protected/seasonally-closed sites (see ZONES + SEASONAL_SPECIES_CLOSURES).
SPECIES = {
    "mutton":      {"label": "Mutton snapper",        "months": [5,6,7,8],     "lunar": "full", "temp": (74,82), "depth": (60,300),
                    "note": "Aggregates on full moons late spring–summer. The main aggregation (Riley's Hump) sits in a YEAR-ROUND no-take reserve."},
    "graysnapper": {"label": "Gray (mangrove) snapper","months": [6,7,8,9],    "lunar": "full", "temp": (78,86), "depth": (30,180),
                    "note": "Summer full/new-moon spawns on reefs and wrecks."},
    "tarpon":      {"label": "Tarpon",                "months": [4,5,6,7],     "lunar": "both", "temp": (74,84), "depth": (10,120),
                    "note": "Catch-and-release fishery in FL. Spawns offshore around new & full moons."},
    "permit":      {"label": "Permit",                "months": [4,5,6,7],     "lunar": "full", "temp": (75,84), "depth": (40,220),
                    "note": "Western Dry Rocks is the key Lower-Keys spawning site and is CLOSED Apr 1–Jul 31."},
    "grouper":     {"label": "Black / gag grouper",   "months": [12,1,2,3],    "lunar": "full", "temp": (68,78), "depth": (60,330),
                    "note": "Atlantic shallow-water grouper is CLOSED Jan 1–Apr 30. Gag has extra restrictions."},
}

# ===================== CONSERVATION / REGULATIONS =====================
# Approximate centres + radii (nautical miles) for no-take / seasonal sanctuary zones
# relevant to the Keys. Coordinates are APPROXIMATE — confirm exact boundaries on
# official charts / the Fish Rules app. Sources: NOAA Florida Keys NMS, FWC (2026).
ZONES = [
    {"name": "Western Sambo Ecological Reserve",      "lat": 24.480, "lon": -81.720, "r": 1.7, "kind": "no_take"},
    {"name": "Western Dry Rocks (spawning closure)",  "lat": 24.448, "lon": -81.926, "r": 0.7, "kind": "seasonal", "months": [4,5,6,7]},
    {"name": "Tortugas South / Riley's Hump Reserve", "lat": 24.490, "lon": -83.100, "r": 4.5, "kind": "no_take"},
    {"name": "Tortugas North Ecological Reserve",     "lat": 24.660, "lon": -83.100, "r": 5.4, "kind": "no_take"},
    {"name": "Dry Tortugas NP Research Natural Area", "lat": 24.630, "lon": -82.870, "r": 3.6, "kind": "no_take"},
]
# FWC/NOAA seasonal HARVEST closures by species group (Atlantic state waters; 2026).
# VERIFY every season in the Fish Rules app — dates shift year to year.
SEASONAL_SPECIES_CLOSURES = {
    "grouper": {"months": [1,2,3,4],   "what": "Atlantic shallow-water grouper closed Jan 1–Apr 30 (spawning)"},
    "mutton":  {"months": [],          "what": "Check current mutton snapper size/bag + any spawning-season rules"},
    "permit":  {"months": [],          "what": "Permit harvest closed in the Keys during the spring spawn; special-permit zone rules"},
}

REG_LINKS = {
    "FWC saltwater regs": "https://myfwc.com/fishing/saltwater/recreational/",
    "Fish Rules app":     "https://www.fishrulesapp.com/",
    "FL Keys NMS zones":  "https://floridakeys.noaa.gov/zones/",
}

# ============================ HELPERS ============================
def moon_fraction(dt):
    synodic = 29.53058867
    known_new = datetime.datetime(2000, 1, 6, 18, 14, tzinfo=datetime.timezone.utc)
    days = (dt - known_new).total_seconds() / 86400.0
    p = (days % synodic) / synodic
    return p + 1 if p < 0 else p

def lunar_match(p, kind):
    if kind == "full": return abs(p - 0.5) < 0.12
    if kind == "new":  return p < 0.12 or p > 0.88
    if kind == "both": return (p < 0.12 or p > 0.88) or abs(p - 0.5) < 0.12
    return False

def haversine_nm(la1, lo1, la2, lo2):
    R = 3440.065  # nautical miles
    p = math.pi / 180
    dla = (la2 - la1) * p; dlo = (lo2 - lo1) * p
    a = math.sin(dla/2)**2 + math.cos(la1*p)*math.cos(la2*p)*math.sin(dlo/2)**2
    return 2 * R * math.asin(math.sqrt(a))

def zone_flags(lat, lon, month):
    """Return conservation flags for a point: protected / closed-now / seasonal-open."""
    flags = []
    for z in ZONES:
        if haversine_nm(lat, lon, z["lat"], z["lon"]) <= z["r"]:
            if z["kind"] == "no_take":
                flags.append(("PROTECTED", f"No-take reserve: {z['name']} — do not fish"))
            elif z["kind"] == "seasonal":
                if month in z.get("months", []):
                    flags.append(("CLOSED", f"Seasonal closure active: {z['name']}"))
                else:
                    flags.append(("NOTE", f"Seasonal zone (open now): {z['name']}"))
    return flags

def species_closed(species, month):
    c = SEASONAL_SPECIES_CLOSURES.get(species)
    if c and month in c.get("months", []):
        return c["what"]
    return None

# ============================ DATA ============================
_DIMS_CACHE = {}
def dataset_dims(server, dataset, var, timeout=30):
    """Read the dataset's structure (.dds) once and return var's dimension names in order,
    e.g. ['time','altitude','latitude','longitude']. None if it can't be read.
    Needed because ERDDAP rejects a request whose bracket count doesn't match: the VIIRS
    chlorophyll products carry a dummy 'altitude' axis that the old MODIS ones didn't."""
    key = (server, dataset, var)
    if key not in _DIMS_CACHE:
        dims = None
        try:
            r = requests.get(f"{server}/griddap/{dataset}.dds", timeout=timeout)
            if r.status_code == 200:
                m = re.search(rf"\b{re.escape(var)}((?:\[\w+ = \d+\])+)", r.text)
                if m:
                    dims = re.findall(r"\[(\w+) = \d+\]", m.group(1))
            else:
                print(f"  ! {dataset}.dds: HTTP {r.status_code}")
        except Exception as e:
            print(f"  ! {dataset}.dds: {e}")
        _DIMS_CACHE[key] = dims
    return _DIMS_CACHE[key]

def build_selector(var, dims, box, stride, time_sel):
    """One bracket per dimension, in the dataset's own order. Extra axes (altitude, depth,
    level…) take index 0. Index and value forms may be mixed per dimension in ERDDAP."""
    b = box
    parts = []
    for d in dims:
        dl = d.lower()
        if dl == "time":
            parts.append(time_sel or "[last]")
        elif dl in ("latitude", "lat"):
            parts.append(f"[({b['latmax']}):{stride}:({b['latmin']})]")
        elif dl in ("longitude", "lon"):
            parts.append(f"[({b['lonmin']}):{stride}:({b['lonmax']})]")
        else:
            parts.append("[0]")
    return var + "".join(parts)

def fetch_grid(dataset, var, box, stride=1, timeout=90, server=ERDDAP_PFEG, time_sel="[(last)]"):
    """Pull one variable over the box. Returns {"grid","lats","lons","time","n_times"} or None.
    If time_sel spans several slices (e.g. "[last-6:last]") the per-pixel MEDIAN over time
    is returned — a cheap cloud-hole filler. "time" is the newest slice's ISO date."""
    dims = dataset_dims(server, dataset, var)
    if not dims:   # structure unreadable — fall back to the classic layout
        dims = (["time"] if time_sel else []) + ["latitude", "longitude"]
    sel = build_selector(var, dims, box, stride, time_sel)
    try:
        r = None
        for attempt in (sel, sel.replace("[0]", "[(0.0)]")):   # 2nd form only if the server dislikes mixed index/value brackets
            url = f"{server}/griddap/{dataset}.json?{urllib.parse.quote(attempt, safe='()[]:.,-')}"
            r = requests.get(url, timeout=timeout)
            if r.status_code == 200 or "[0]" not in sel or r.status_code != 400:
                break
        if r.status_code != 200:
            why = re.sub(r"\s+", " ", r.text)[:160].strip()   # ERDDAP says WHY in the body — keep it in the log
            print(f"  ! {dataset}: HTTP {r.status_code} {why}")
            return None
        t = r.json().get("table")
        if not t:
            return None
        cols, rows = t["columnNames"], t["rows"]
        ci, cj, cv = cols.index("latitude"), cols.index("longitude"), cols.index(var)
        ct = cols.index("time") if "time" in cols else None
        lats = sorted({row[ci] for row in rows}, reverse=True)
        lons = sorted({row[cj] for row in rows})
        li = {v: i for i, v in enumerate(lats)}
        lj = {v: i for i, v in enumerate(lons)}
        times = sorted({row[ct] for row in rows}) if ct is not None else [None]
        ti = {v: i for i, v in enumerate(times)}
        stack = np.full((len(times), len(lats), len(lons)), np.nan)
        for row in rows:
            if row[cv] is not None:
                stack[ti[row[ct]] if ct is not None else 0, li[row[ci]], lj[row[cj]]] = row[cv]
        with warnings.catch_warnings():                 # all-NaN pixels are expected (clouds/land)
            warnings.simplefilter("ignore", category=RuntimeWarning)
            g = np.nanmedian(stack, axis=0) if len(times) > 1 else stack[0]
        return {"grid": g, "lats": np.array(lats), "lons": np.array(lons),
                "time": times[-1], "n_times": len(times)}
    except Exception as e:
        print(f"  ! {dataset}: {e}")
        return None

def age_days(iso):
    """Age of an ERDDAP ISO timestamp in days (float), or None if unparseable."""
    if not iso:
        return None
    try:
        t = datetime.datetime.fromisoformat(str(iso).replace("Z", "+00:00"))
        return (datetime.datetime.now(datetime.timezone.utc) - t).total_seconds() / 86400.0
    except ValueError:
        return None

def load_from_ladder(kind, sources, box, stride):
    """Walk a SST/CHL ladder; return the first FRESH grid (with its source attached) or None.
    Every rung's outcome is printed so the run log shows exactly what fed the score."""
    for s in sources:
        tsel = f"[last-{s['days'] - 1}:last]" if s["days"] > 1 else "[(last)]"
        g = fetch_grid(s["dataset"], s["var"], box, stride, server=s["server"], time_sel=tsel)
        if not g:
            print(f"  {kind}: {s['label']} — no answer, next rung")
            continue
        age = age_days(g["time"])
        if age is None:
            print(f"  {kind}: {s['label']} — image has no date, next rung")
            continue
        if age > s["max_age_days"]:
            print(f"  {kind}: {s['label']} — newest image {g['time'][:10]} is {age:.0f} days old "
                  f"(limit {s['max_age_days']}) — STALE, next rung")
            continue
        comp = f", {g['n_times']}-day composite" if g["n_times"] > 1 else ""
        print(f"  {kind}: {s['label']} — newest image {g['time'][:10]} ({age:.1f} days old{comp})  OK")
        g["source"] = s
        g["date"] = g["time"][:10]
        return g
    print(f"  !! {kind}: every source in the ladder was down or stale")
    return None

def resample_to(src, tgt_lats, tgt_lons):
    """Nearest-neighbour resample src grid onto target lat/lon axes."""
    g, sl, so = src["grid"], src["lats"], src["lons"]
    out = np.full((len(tgt_lats), len(tgt_lons)), np.nan)
    ii = [int(np.argmin(np.abs(sl - la))) for la in tgt_lats]
    jj = [int(np.argmin(np.abs(so - lo))) for lo in tgt_lons]
    for a, i in enumerate(ii):
        for b, j in enumerate(jj):
            out[a, b] = g[i, j]
    return out

# ----- structure (known wrecks / artificial reefs) -----
NAME_KEYS = ["Name", "NAME", "name", "Reef_Name", "REEF_NAME", "AR_NAME", "SITE_NAME",
             "Site_Name", "FEATURE_NAME", "vesselterm", "VESSELTERM", "feature_type",
             "FEATURE_TYPE", "sondat", "COUNTY", "County"]

def _pick_name(props, default):
    for k in NAME_KEYS:
        v = props.get(k)
        if isinstance(v, str) and v.strip() and v.strip().lower() not in ("null", "unknown", "n/a"):
            return v.strip()[:48]
    return default

def fetch_arcgis_points(layer_url, box, kind, timeout=45, cap=2500):
    """Query an ArcGIS REST layer with the Keys envelope, return [{name,lat,lon,kind}]."""
    b = box
    params = {
        "where": "1=1",
        "geometry": f"{b['lonmin']},{b['latmin']},{b['lonmax']},{b['latmax']}",
        "geometryType": "esriGeometryEnvelope", "inSR": "4326",
        "spatialRel": "esriSpatialRelIntersects", "outFields": "*",
        "returnGeometry": "true", "outSR": "4326", "f": "geojson",
    }
    url = f"{layer_url}/query?{urllib.parse.urlencode(params)}"
    try:
        r = requests.get(url, timeout=timeout)
        if r.status_code != 200:
            print(f"  ! {kind}: HTTP {r.status_code}")
            return []
        feats = r.json().get("features", []) or []
        out = []
        for f in feats[:cap]:
            geom = (f.get("geometry") or {}).get("coordinates")
            if not geom or len(geom) < 2 or geom[0] is None or geom[1] is None:
                continue
            lon, lat = float(geom[0]), float(geom[1])
            out.append({"name": _pick_name(f.get("properties", {}) or {}, kind.title()),
                        "lat": round(lat, 5), "lon": round(lon, 5), "kind": kind})
        print(f"  {kind}: {len(out)} points in box")
        return out
    except Exception as e:
        print(f"  ! {kind}: {e}")
        return []

def gather_structure(box):
    pts = []
    if CONFIG["use_structure"]:
        pts += fetch_arcgis_points(FWC_REEFS, box, "reef")
        pts += fetch_arcgis_points(NOAA_WRECKS, box, "wreck")
    return pts

def _hav_nm_arr(lat0, lon0, LAT, LON):
    R = 3440.065
    p = math.pi / 180
    dla = (LAT - lat0) * p; dlo = (LON - lon0) * p
    a = np.sin(dla/2)**2 + math.cos(lat0*p) * np.cos(LAT*p) * np.sin(dlo/2)**2
    return 2 * R * np.arcsin(np.sqrt(np.clip(a, 0, 1)))

def structure_boost(points, lats, lons, radius_nm):
    """0..1 grid: cells near a wreck/reef get a smooth bump (max over all points)."""
    g = np.zeros((len(lats), len(lons)))
    if not points:
        return g
    LON, LAT = np.meshgrid(np.asarray(lons, float), np.asarray(lats, float))
    for s in points:
        d = _hav_nm_arr(s["lat"], s["lon"], LAT, LON)
        bump = np.where(d <= radius_nm * 2.5, np.exp(-(d / radius_nm) ** 2), 0.0)
        np.maximum(g, bump, out=g)
    return g

def nearest_structure(lat, lon, points):
    best, bd = None, 1e9
    for s in points:
        d = haversine_nm(lat, lon, s["lat"], s["lon"])
        if d < bd:
            bd, best = d, s
    return best, bd

# ============================ MODEL ============================
def grad_mag(a):
    filled = np.where(np.isnan(a), np.nanmean(a), a)
    gy, gx = np.gradient(filled)
    m = np.hypot(gx, gy)
    m[np.isnan(a)] = np.nan
    return m

def norm98(a):
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return np.zeros_like(a)
    hi = np.percentile(finite, 98) or 1.0
    return np.clip(np.nan_to_num(a) / hi, 0, 1)

def current_fields(U, V, lats, lons):
    """From eastward/northward velocity grids compute, on the native grid:
       speed (m/s), surface CONVERGENCE (>0 = water piling up: bait & weed gather),
       and an EDGE field (gradient of speed: the 'wall' of a strong current).
    Returns (speed, convergence, edge) as same-shape arrays, NaN where no data."""
    nan = np.isnan(U) | np.isnan(V)
    speed = np.hypot(U, V)
    # physical spacing in metres (lat axis may be descending — np.gradient handles signed coords)
    mean_lat = float(np.nanmean(lats))
    y_m = np.asarray(lats, float) * 110540.0
    x_m = np.asarray(lons, float) * 111320.0 * math.cos(math.radians(mean_lat))
    Uf = np.where(nan, np.nanmean(U), U)
    Vf = np.where(nan, np.nanmean(V), V)
    if len(lats) < 2 or len(lons) < 2:
        conv = np.zeros_like(speed)
    else:
        dUdy, dUdx = np.gradient(Uf, y_m, x_m)
        dVdy, dVdx = np.gradient(Vf, y_m, x_m)
        conv = -(dUdx + dVdy)            # convergence = -divergence
    edge = grad_mag(speed)
    speed[nan] = np.nan; conv[nan] = np.nan; edge[nan] = np.nan
    return speed, conv, edge

def weed_anomaly(afai):
    """Sargassum signal from AFAI: how far above the local background each cell sits.
    AFAI saturates under sun glint / cloud, so treat as best-effort, not gospel."""
    med = np.nanmedian(afai)
    a = afai - (med if np.isfinite(med) else 0.0)
    return np.clip(a, 0, None)

def analyze():
    box, stride = CONFIG["box"], CONFIG["stride"]
    print("Fetching satellite + bathymetry grids over the Florida Keys…")
    chl = load_from_ladder("chl", CHL_SOURCES, box, stride)
    sst = load_from_ladder("SST", SST_SOURCES, box, stride)
    bathy = fetch_grid(CONFIG["bathy_dataset"], CONFIG["bathy_var"], box,
                       max(CONFIG["bathy_stride"], 1), time_sel="")
    if not bathy and CONFIG["bathy_dataset"] != "etopo180":
        print("  bathy: srtm15plus unavailable — falling back to ETOPO")
        bathy = fetch_grid("etopo180", "altitude", box, max(stride, 1), time_sel="")

    if not chl or not sst:
        missing = " and ".join(k for k, v in (("chlorophyll", chl), ("sea temperature", sst)) if not v)
        msg = (f"BiteCast: no FRESH satellite {missing} today — every source in the ladder was "
               "down or older than its limit. No edge report this run (the run log lists each "
               "source and its image date).")
        print(msg)
        return push(msg, title="BiteCast — data unavailable", tags="warning,fish")
    data_line = (f"Data: SST {sst['source']['label']} ({sst['date']}) · "
                 f"chl {chl['source']['label']} ({chl['date']})")

    lats, lons = chl["lats"], chl["lons"]
    C = chl["grid"]
    S = resample_to(sst, lats, lons)
    D = resample_to(bathy, lats, lons) if bathy else None  # metres, negative below sea level

    # --- fronts (coincident temp + colour breaks) ---
    sstF = norm98(grad_mag(S))
    chlF = norm98(grad_mag(np.log(np.clip(C, 1e-3, None))))
    front = np.sqrt(sstF * chlF)

    # --- structure + species depth band from bathymetry ---
    sp = SPECIES[CONFIG["target_species"]]
    if D is not None:
        slope = norm98(grad_mag(D))                 # drop-offs, ledges, humps
        depth_ft = np.abs(D) * 3.28084
        lo_d, hi_d = sp["depth"]
        band = np.where((depth_ft >= lo_d) & (depth_ft <= hi_d), 1.0,
                        np.clip(1 - np.minimum(np.abs(depth_ft - lo_d), np.abs(depth_ft - hi_d)) / 150.0, 0, 1))
        on_water = D < 0                            # mask land
    else:
        slope = np.zeros_like(front); band = np.ones_like(front); on_water = np.ones_like(front, bool)

    # --- ocean currents: Florida-Current edge + convergence (bait & weed accumulate here) ---
    cur_kt = convN = edgeN = None
    if CONFIG["use_currents"]:
        U = fetch_grid(CONFIG["current_dataset"], "u_current", box, 1, server=ERDDAP_CW)
        V = fetch_grid(CONFIG["current_dataset"], "v_current", box, 1, server=ERDDAP_CW)
        if U and V:
            speed_ms, conv, edge = current_fields(U["grid"], V["grid"], U["lats"], U["lons"])
            cur_kt = resample_to({"grid": speed_ms, "lats": U["lats"], "lons": U["lons"]}, lats, lons) * 1.943844
            convN  = norm98(resample_to({"grid": conv, "lats": U["lats"], "lons": U["lons"]}, lats, lons))
            edgeN  = norm98(resample_to({"grid": edge, "lats": U["lats"], "lons": U["lons"]}, lats, lons))
            print(f"  currents: up to {np.nanmax(cur_kt):.1f} kt; edge + convergence layered in")
        else:
            print("  ! currents unavailable this run — score falls back to the other layers")

    # --- sargassum / weed lines (AFAI floating-algae index) ---
    weedN = afai = None
    if CONFIG["use_sargassum"]:
        A = fetch_grid(CONFIG["sargassum_dataset"], "AFAI", box, 1, server=ERDDAP_AOML)
        if A:
            afai  = resample_to({"grid": A["grid"], "lats": A["lats"], "lons": A["lons"]}, lats, lons)
            weedN = norm98(resample_to({"grid": weed_anomaly(A["grid"]), "lats": A["lats"], "lons": A["lons"]}, lats, lons))
            print("  sargassum (AFAI) layered in")
        else:
            print("  ! sargassum (AFAI) unavailable this run — skipping the weed signal")

    # --- structure: known wrecks + artificial reefs (exact coordinates) ---
    structures = gather_structure(box)
    reefN = structure_boost(structures, lats, lons, CONFIG["structure_radius_nm"]) if structures else None

    # --- modular fusion: weight only the layers that loaded, then renormalise ---
    layers = [(CONFIG["w_front"], front)]
    if D is not None:
        layers += [(CONFIG["w_structure"], slope), (CONFIG["w_depth"], band)]
    if reefN is not None:
        layers += [(CONFIG["w_reef"], reefN)]
    if edgeN is not None:
        layers += [(CONFIG["w_current"], edgeN), (CONFIG["w_converge"], convN)]
    if weedN is not None:
        layers += [(CONFIG["w_weed"], weedN)]
    tot_w = sum(w for w, _ in layers) or 1.0
    fish = sum(w * arr for w, arr in layers) / tot_w
    fish = np.where(on_water, fish, np.nan)

    # --- spawning suitability (season + moon + temp + depth band) ---
    now = datetime.datetime.now(datetime.timezone.utc)
    month = now.month
    p = moon_fraction(now)
    in_season = month in sp["months"]
    moon_ok = lunar_match(p, sp["lunar"])
    sst_f = np.where(np.isnan(S), np.nan, S * 9/5 + 32)
    lo_t, hi_t = sp["temp"]
    temp_ok = (sst_f >= lo_t) & (sst_f <= hi_t)
    spawn = (front * band * (temp_ok.astype(float))) if (in_season and moon_ok) else np.zeros_like(front)

    # --- rank hotspots, attach conservation flags ---
    flat = [(fish[i, j], i, j) for i in range(len(lats)) for j in range(len(lons))
            if np.isfinite(fish[i, j]) and fish[i, j] >= CONFIG["min_edge_score"]]
    flat.sort(reverse=True)

    spots = []
    for score, i, j in flat[:80]:
        lat, lon = float(lats[i]), float(lons[j])
        fl = zone_flags(lat, lon, month)
        protected = any(k == "PROTECTED" or k == "CLOSED" for k, _ in fl)
        cur = None if cur_kt is None or np.isnan(cur_kt[i, j]) else round(float(cur_kt[i, j]), 2)
        cvg = None if convN  is None or np.isnan(convN[i, j])  else round(float(convN[i, j]), 3)
        wd  = None if weedN  is None or np.isnan(weedN[i, j])  else round(float(weedN[i, j]), 3)
        tags = []
        if edgeN is not None and edgeN[i, j] >= 0.5: tags.append("current edge")
        if cvg is not None and cvg >= 0.5:           tags.append("convergence/accumulation")
        if wd  is not None and wd  >= 0.4:           tags.append("likely weed line")
        near = nd = None
        if structures:
            ns, dist = nearest_structure(lat, lon, structures)
            if ns and dist <= 3.0:
                near, nd = f"{ns['name']} ({ns['kind']})", round(dist, 2)
                if dist <= CONFIG["structure_radius_nm"]:
                    tags.append("on structure")
        spots.append({
            "lat": round(lat, 4), "lon": round(lon, 4),
            "fish_score": round(float(score), 3),
            "front": round(float(front[i, j]), 3),
            "structure": round(float(slope[i, j]), 3),
            "sst_f": None if np.isnan(sst_f[i, j]) else round(float(sst_f[i, j]), 1),
            "chl": None if np.isnan(C[i, j]) else round(float(C[i, j]), 3),
            "depth_ft": None if D is None or np.isnan(D[i, j]) else round(float(abs(D[i, j]) * 3.28084)),
            "cur_kt": cur, "convergence": cvg, "weed": wd, "tags": tags,
            "near_structure": near, "near_structure_nm": nd,
            "spawn_score": round(float(spawn[i, j]), 3),
            "dist_home_nm": round(haversine_nm(CONFIG["home"]["lat"], CONFIG["home"]["lon"], lat, lon), 1),
            "flags": [f"{k}: {m}" for k, m in fl],
            "protected": protected,
        })

    # a one-line read on the current/weed regime for the push
    regime = [data_line]
    if cur_kt is not None:
        regime.append(f"Currents live — up to {np.nanmax(cur_kt):.1f} kt in the box "
                      f"(Florida-Current edge & eddies factored in).")
    if weedN is not None:
        regime.append("Sargassum/weed signal present in cloud-free pixels."
                      if int(np.sum(weedN >= 0.4)) else "Little sargassum in the latest clear pixels.")
    if structures:
        nr = sum(1 for s in structures if s["kind"] == "reef")
        nw = sum(1 for s in structures if s["kind"] == "wreck")
        regime.append(f"Structure: {nr} reefs + {nw} wrecks in range (display only — verify, not for nav).")

    # keep top, but never recommend protected/closed water — flag it instead
    openable = [s for s in spots if not s["protected"]][:CONFIG["top_n"]]
    blocked = [s for s in spots if s["protected"]][:4]

    meta = {
        "generated_at": now.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sst": {"source": sst["source"]["label"], "dataset": sst["source"]["dataset"], "image_date": sst["date"]},
        "chl": {"source": chl["source"]["label"], "dataset": chl["source"]["dataset"], "image_date": chl["date"],
                "composite_days": chl["n_times"]},
        "grid_cells": int(np.isfinite(fish).sum()),
        "layers": [name for name, on in (("front", True), ("bathymetry", D is not None), ("reefs_wrecks", reefN is not None),
                                         ("currents", edgeN is not None), ("sargassum", weedN is not None)) if on],
        "target_species": CONFIG["target_species"],
    }
    write_geojson(spots[:40], meta)
    write_structure_geojson(structures)
    return report_and_push(openable, blocked, sp, in_season, moon_ok, p, month, regime)

# ============================ OUTPUT ============================
def write_geojson(spots, meta=None):
    # "meta" stamps the file with when it was made and which images fed it, so a stale
    # file can never pass for a fresh one (the web app can show these dates).
    fc = {"type": "FeatureCollection", "meta": meta or {}, "features": []}
    for s in spots:
        fc["features"].append({
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [s["lon"], s["lat"]]},
            "properties": s,
        })
    with open("hotspots.geojson", "w") as f:
        json.dump(fc, f, indent=1)
    print(f"Wrote hotspots.geojson ({len(spots)} features)")

def write_structure_geojson(points):
    fc = {"type": "FeatureCollection",
          "note": "FWC artificial reefs + NOAA AWOIS wrecks. Display/planning only — NOT for navigation.",
          "features": []}
    for s in points:
        fc["features"].append({
            "type": "Feature",
            "geometry": {"type": "Point", "coordinates": [s["lon"], s["lat"]]},
            "properties": {"name": s["name"], "kind": s["kind"]},
        })
    with open("structure.geojson", "w") as f:
        json.dump(fc, f, indent=1)
    print(f"Wrote structure.geojson ({len(points)} reefs/wrecks)")

NTFY_URL = os.environ.get("NTFY_URL", "https://ntfy.sh")   # override for a self-hosted ntfy
NTFY_MAX_BYTES = 3900   # ntfy turns bodies over 4,096 bytes into a .txt attachment; stay under it
_PRIORITY = {"min": 1, "low": 2, "default": 3, "high": 4, "max": 5, "urgent": 5}
TOPIC_RE = re.compile(r"^[-_A-Za-z0-9]{1,64}$")   # ntfy's own rule for topic names

def _clean_topic(raw):
    """Accept the topic as typed, or as pasted from the app: 'ntfy.sh/topic',
    'https://ntfy.sh/topic', with stray spaces or a trailing slash."""
    t = (raw or "").strip()
    if "/" in t:
        t = t.rstrip("/").rsplit("/", 1)[-1].strip()
    return t

def _fit_ntfy(body):
    """Trim the body at a line boundary so the phone shows text, not an attachment."""
    if len(body.encode("utf-8")) <= NTFY_MAX_BYTES:
        return body
    out, used = [], 0
    for line in body.split("\n"):
        n = len(line.encode("utf-8")) + 1
        if used + n > NTFY_MAX_BYTES - 24:
            break
        out.append(line); used += n
    return "\n".join(out) + "\n…(trimmed to fit)"

def push(body, title="BiteCast", tags="fish", priority="default"):
    """Send one ntfy notification. Returns True only when ntfy accepted it.

    Uses ntfy's JSON publish endpoint: the title travels in the JSON body, so
    characters like '—' and '•' are fine. (The old header-based call died on the
    em dash in the title before the request ever left the machine, and the error
    was swallowed, so the run looked green while the phone got nothing.)"""
    topic = _clean_topic(CONFIG["ntfy_topic"])
    if not topic or "CHANGE-ME" in topic:
        print("!! NTFY_TOPIC is not set (repo Settings → Secrets → Actions). Nothing was pushed.\n"
              "   Report that would have been sent:\n" + body)
        return False
    if not TOPIC_RE.match(topic):
        bad = sorted({c for c in topic if not re.match(r"[-_A-Za-z0-9]", c)})
        print(f"!! NTFY_TOPIC is not a valid ntfy topic ({len(topic)} chars, contains {bad}). "
              "Allowed: letters, digits, '-' and '_', up to 64 chars. Put ONLY the topic word in the "
              "secret, not the URL. Nothing was pushed.")
        return False
    payload = {
        "topic": topic,
        "title": title,
        "message": _fit_ntfy(body),
        "tags": [t.strip() for t in tags.split(",") if t.strip()],
        "priority": _PRIORITY.get(str(priority).lower(), 3),
    }
    try:
        r = requests.post(NTFY_URL.rstrip("/") + "/", json=payload, timeout=30)
    except Exception as e:
        print(f"!! ntfy push FAILED before delivery: {type(e).__name__}: {e}")
        return False
    if r.status_code != 200:
        print(f"!! ntfy REJECTED the push: HTTP {r.status_code} {r.text[:200].strip()}")
        return False
    try:
        mid = r.json().get("id", "?")
    except ValueError:
        mid = "?"
    print(f"Pushed to ntfy (message id {mid}).")   # topic deliberately not echoed into the public log
    return True

def ping():
    """Prove the phone channel end to end without touching any data source.
    Run:  python analyze.py --ping   (or the 'ping' choice on the Actions tab)."""
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    ok = push(f"BiteCast is wired to this phone.\nSent {now} from the cloud analyzer.\n"
              f"Home: {CONFIG['home']['name']} · target: {SPECIES[CONFIG['target_species']]['label']}",
              title="BiteCast — test ping", tags="white_check_mark,fish", priority="default")
    print("PING OK — check your phone." if ok else "PING FAILED — see the message above.")
    return ok

def report_and_push(openable, blocked, sp, in_season, moon_ok, p, month, regime=None):
    moon_pct = round((1 - math.cos(2*math.pi*p)) / 2 * 100)
    lines = []
    if regime:
        lines += regime + [""]
    if openable:
        lines.append("Top fishable edges (ranked by strength):")
        for s in sorted(openable, key=lambda x: -x["fish_score"]):
            mark = ("  [" + ", ".join(s["tags"]) + "]") if s.get("tags") else ""
            struct = (f"  ~{s['near_structure_nm']}nm from {s['near_structure']}"
                      if s.get("near_structure") else "")
            lines.append(f"• {s['lat']}N {abs(s['lon'])}W  edge {int(s['fish_score']*100)}/100"
                         f"{', '+str(s['sst_f'])+'F' if s['sst_f'] else ''}"
                         f"{', chl '+str(s['chl']) if s['chl'] else ''}"
                         f"{', '+str(s['cur_kt'])+'kt' if s.get('cur_kt') else ''}"
                         f"{', '+str(s['depth_ft'])+'ft' if s['depth_ft'] else ''}"
                         f"  (~{s['dist_home_nm']} nm from {CONFIG['home']['name']})" + struct + mark)
    else:
        lines.append("No strong open-water edges in the latest cloud-free data this run.")

    # spawning + conservation block
    status = "PEAK window" if (in_season and moon_ok) else ("IN SEASON" if in_season else "off-season")
    lines.append("")
    lines.append(f"{sp['label']} spawning: {status}  (moon {moon_pct}% illum).")
    closed = species_closed(CONFIG["target_species"], month)
    if closed:
        lines.append(f"⚠ HARVEST CLOSED: {closed}. Catch-and-release / observe only — do not keep.")
    lines.append(f"Note: {sp['note']}")
    if blocked:
        lines.append("")
        lines.append("Strong edges that fall in PROTECTED/CLOSED water — do NOT fish these:")
        for s in blocked:
            lines.append(f"• {s['lat']}N {abs(s['lon'])}W — {' / '.join(s['flags'])}")
    lines.append("")
    lines.append("Confirm all rules in the Fish Rules app before fishing. Respect every closure.")

    body = "\n".join(lines)
    # alert priority high only when there's genuinely good open water and nothing's closed
    pr = "high" if (openable and not closed) else "default"
    print("\n" + body + "\n")
    return push(body, title="BiteCast — Keys edge report", priority=pr)

if __name__ == "__main__":
    # A green Actions run now means "the phone got the report". If the push fails
    # (no topic, bad reply, network), the run exits 2 and shows red on the Actions tab;
    # the geojson files are still written and committed by the workflow's next step.
    mode = "ping" if ("--ping" in sys.argv[1:] or os.environ.get("BITECAST_MODE", "").lower() == "ping") else "scan"
    delivered = ping() if mode == "ping" else analyze()
    if not delivered:
        print("\nRESULT: phone NOT notified — failing this run so it shows red on the Actions tab.")
        sys.exit(2)
    print("\nRESULT: phone notified.")
