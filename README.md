# BiteCast — Florida Keys cloud analyzer

A free, scheduled "edge finder + spawning + conservation" engine for the Keys.
It runs in the cloud (GitHub Actions) a few times a day — **your PC can be off** —
pulls free satellite + bathymetry data, finds the productive water, and texts your
phone, while **refusing to point you at closed or protected water**.

## What it does each run
1. Pulls **SST**, **chlorophyll**, and **bathymetry** grids over the Keys (NOAA ERDDAP — free, no key). SST and chlorophyll each come from a **ladder** of datasets (`SST_SOURCES` / `CHL_SOURCES` in `analyze.py`): the run tries them in order of quality, reads the date of the newest image, **refuses anything older than that source's limit**, and logs which one it used. Daily 4 km chlorophyll images are stacked into a **7-day median composite** to fill cloud holes. (The MODIS-Aqua feeds this tool started on stopped updating in 2022 — a fixed dataset id is a trap; a ladder with a freshness check isn't.)
2. Detects **temperature fronts** and **colour breaks**, and the places where they **coincide** (the high-odds water).
3. Adds **bottom structure** (drop-offs / ledges / humps from the depth gradient) and a **species depth band** — now off **~500 m bathymetry** (SRTM15+, which folds in NOAA's Coastal Relief Model near the coast) instead of the old ~1.85 km grid.
4. Pulls **known wrecks & artificial reefs** (FWC's 4,400+ reef points + NOAA AWOIS wrecks) and (a) boosts hotspots that sit **on structure**, (b) tags every top spot with the **nearest wreck/reef and its distance**, and (c) writes **`structure.geojson`** so you can drop the reef/wreck pins on a map or chartplotter.
5. Pulls **ocean currents** (near-real-time altimetry geostrophic, NOAA CoastWatch) and derives the **current edge** — the "wall" of the Florida Current where pelagics stack — plus **surface convergence**, where water (and bait/weed) piles up.
6. Pulls **sargassum** via the USF/AOML **AFAI** floating-algae index and flags **likely weed lines** (mahi / tripletail / bait magnets).
7. **Fuses** every layer that loaded into one *fish-suitability* score (the weighting **auto-renormalises**, so if any source is briefly unavailable the score still works).
8. Layers a **spawning-window** model for your target species (season + moon phase + temperature + depth).
9. Applies a **conservation filter**: any hotspot inside a no-take reserve or an active seasonal closure is **flagged and withheld from the "go fish" list**, with the reason.
10. Pushes a summary to your phone via **ntfy** (each spot tagged *current edge / convergence / weed line / on structure*, with the nearest reef/wreck), and writes three files the web app reads from the same GitHub Pages folder: **`hotspots.geojson`** (the fused hotspots + protected-zone flags, with a `meta` block), **`structure.geojson`** (reefs/wrecks) and **`satgrid.json`** (chlorophyll + sea temperature on the analysis grid, for the app's condition cards). Browsers can't read NOAA's ERDDAP directly (no CORS header), so the app never talks to NOAA itself — it shows what the last scan found and says when that was.

### Data sources for the new signals (all free, no key)
- **Bathymetry:** `srtm15plus` on `coastwatch.pfeg.noaa.gov/erddap` (var `z`, ~500 m). For the **Atlantic side at ~90 m**, swap `bathy_dataset` to a NOAA Coastal Relief Model volume (e.g. `usgsCeCrm2`, Atlantic Southeast). Note: the *fused* grid runs at the finest satellite layer available (~2 km when MUR 1 km SST is in play, thinned ×2); the coarser layer is interpolated smoothly onto it. Gradients are computed only where every neighbour is real water, so coastlines and cloud edges no longer pose as fronts, and hotspots closer than `min_separation_nm` (3 nm) to a stronger one are dropped so the list shows distinct edges.
- **Structure:** FWC `Artificial_Reef_Locations_in_Florida/MapServer/12` on `gis.myfwc.com` + NOAA AWOIS wrecks on `maps.nccs.nasa.gov` (both ArcGIS REST, queried as GeoJSON). **Display/planning only — neither agency verifies every point, and FWC states plainly: do NOT use for navigation.** Tune influence with `w_reef` and `structure_radius_nm`; set `use_structure=False` to skip.
- **Currents:** `noaacwBLENDEDNRTcurrentsDaily` on `coastwatch.noaa.gov/erddap` — blended altimetry geostrophic, near-real-time, 0.25° (`u_current`/`v_current`, m/s). Coarse (~25 km) but nails the Florida-Current edge and large eddies.
- **Sargassum:** `noaa_aoml_atlantic_oceanwatch_AFAI_7D` on `cwcgom.aoml.noaa.gov/erddap` — USF/AOML 7-day AFAI (~1.6 km). Best-effort: AFAI saturates under sun glint and cloud, so treat weed flags as a strong hint, not a guarantee.
- Tune everything in `CONFIG` (`w_*` weights), or set the `use_*` flags to `False`. Currents + weed help most for **pelagics (mahi, sails, tuna)**; structure + depth dominate for bottom species.

## Setup (~10 min)
1. **Phone:** install the **ntfy** app, subscribe to a unique topic (e.g. `bitecast-keys-7h3xq9`).
2. **GitHub:** create a repo and add these files (keep the `.github/workflows/` path).
3. **Repo → Settings → Secrets and variables → Actions → New repository secret:** name `NTFY_TOPIC`, value = your topic.
4. **Actions tab:** enable workflows, then **Run workflow** with mode **`ping`**. Within a few seconds your phone should show *"BiteCast — test ping"*. Nothing else runs in ping mode — it only proves the phone channel. Then run once more with mode **`scan`** for a real report.
   - **Green run = the phone got the message. Red run = it didn't.** The analyzer exits with an error whenever ntfy did not accept the push (secret missing, ntfy rejected it, network down), so a broken channel can never hide behind a green check mark. The run log says exactly why.
   - Locally: `NTFY_TOPIC=your-topic python analyze.py --ping`.
5. Edit the `CONFIG` block in `analyze.py` to set your **box**, **target_species**, and **home** spot. Adding spots/areas later is just editing this file and committing — no rebuild.

## Conservation — read this
Regulations change constantly. **This is a planning aid, not legal advice.** Always confirm
current rules in the **Fish Rules app** and at **myfwc.com** before fishing, and respect every
closure. Built-in awareness (Atlantic state waters, 2026 — *verify before relying on it*):

- **Atlantic shallow-water grouper** (black, gag, scamp, red, yellowfin, rock/red hind, etc.): **closed Jan 1 – Apr 30** for spawning.
- **Western Dry Rocks** (≈10 nm SW of Key West): **seasonal closure Apr 1 – Jul 31**, continuous transit only, no anchoring — a key permit/multispecies spawning site.
- **Tortugas South / Riley's Hump**, **Tortugas North**, **Western Sambo**, **Dry Tortugas NP Research Natural Area**: **year-round no-take** reserves protecting spawning aggregations. Riley's Hump's closure is why Keys mutton snapper rebounded.
- Snapper–grouper share a **10-fish aggregate** bag limit; Biscayne National Park has its own special rules.

The tool deliberately will **not** recommend fishing protected/closed water — fishing spawning
aggregations strips out the big breeders. Fish the *feeding edges*; let the *spawners* spawn.

## Paid services — are they worth it?
The data underneath all of them is largely **free government satellite data** — paid services
add higher-resolution/composite imagery, human analysis, extra layers (altimetry, currents,
salinity), and chartplotter integration. Rough current pricing:

| Service | ~Cost | Adds over this free tool |
|---|---|---|
| **Hilton's Realtime Navigator** | ~$200/yr per region | SST/chl/altimetry/currents/bathymetry, waypoints, chartplotter (RT-NAV/SAT2NAV) |
| **ROFFS** | ~$36–$65 per custom analysis (≈$1,300 for a full season) | Human oceanographer write-up of your exact zone |
| **FishTrack** | ~$79/yr | Easy SST/chl charts + Buoyweather, phone/tablet friendly |
| **RipCharts** | comparable, often cheaper; Android-friendly | Live charts you can track your GPS position on |
| **SiriusXM Fish Mapping** | marine subscription (needs a compatible receiver/plotter) | Fish-mapping + weather pushed straight to the helm |
| **SatFish** | subscription | High-def imagery + strong bathymetric integration |

**Worth it if** you run offshore often, fish tournaments, or want plotter integration / a human
analyst. **Not needed if** this free tool + the in-app NASA overlays already cover you — you're
using the same satellite feeds, just auto-analyzed for edges.

## Caveats
- Satellite grids are **4–9 km** (chlorophyll) and **1–5 km** (SST) depending on which rung of the ladder was fresh that day — the run log and the push's "Data:" line say which. Expect broad (not pinpoint) edges, and thinner coverage after a cloudy week.
- `hotspots.geojson` carries a `meta` block (`generated_at`, the SST/chl source and image dates, layers that loaded) so a stale file can never pass for a fresh one.
- Coordinates for protected zones are **approximate** — confirm exact boundaries officially.
- The spawning model is **general guidance**, not a guarantee, and never a reason to fish a closure.
- **First-run check:** ERDDAP dataset IDs and variable names occasionally change. If a run logs `currents unavailable` or `sargassum unavailable`, paste the dataset's URL into a browser to confirm it's still live and the variable name still matches (the analyzer prints the dataset IDs it uses). Currents live on `coastwatch.noaa.gov/erddap`, sargassum on `cwcgom.aoml.noaa.gov/erddap`. The score keeps working on whatever layers do load.
