# 3D Gas Dispersion Model

**A Python / Google Colab research and screening tool for atmospheric gas dispersion, gas-mixture properties, release calculations, and preliminary flammability/toxic-exposure screening.**

> [!CAUTION]
> This is a research and screening prototype. It is **not** PHAST, ALOHA, DEGADIS, or a validated replacement for robust consequence-analysis software. Do not use its predictions alone to establish emergency exclusion distances, safety-system design limits, regulatory compliance, or operational emergency-response decisions. Independent validation, scenario-specific source-term assessment, and qualified engineering review are required.

## Contents

- [1. Model overview](#1-model-overview)
- [2. Features and model-selection logic](#2-features-and-model-selection-logic)
- [3. Installation and running in Google Colab](#3-installation-and-running-in-google-colab)
- [4. Quick-start configuration](#4-quick-start-configuration)
- [5. Input configuration reference](#5-input-configuration-reference)
- [6. Model equations and methods](#6-model-equations-and-methods)
- [7. Outputs and interpretation](#7-outputs-and-interpretation)
- [8. Chemical database](#8-chemical-database)
- [9. Assumptions, limitations, and validation](#9-assumptions-limitations-and-validation)
- [10. References](#10-references)

---

## 1. Model overview

The 3D Gas Dispersion Model is a Google Colab/Jupyter application for estimating gas concentrations following a specified atmospheric release. It links a user-defined chemical or fixed-composition gas mixture to a release source term, meteorological conditions, an atmospheric-dispersion engine, receptor calculations, and preliminary hazard screening.

The model uses local wind-aligned coordinates:

- **x**: downwind distance, positive in the mean wind direction;
- **y**: crosswind distance;
- **z**: elevation above ground.

The source is located at `(x = 0, y = source_y_m, z = release_height_m)`. Wind direction is entered as a meteorological **wind-from** direction. Dispersion plots remain in local wind-aligned coordinates; an optional CSV export transforms x-y offsets to east/north coordinates.

### 1.1 Calculation workflow

The calculation sequence is presented below as a numbered text workflow. It uses no diagram or external image file.

**01. Define the scenario**

- Enter the release location, release elevation, computational domain, observation time, and receptor locations.
- Set the meteorological conditions: reference wind speed and height, wind direction, atmospheric stability, surface roughness, temperature, and pressure.

**02. Select the chemical mixture**

- Select one or more chemicals from the built-in chemical database.
- Enter each component fraction on a mole or mass basis.
- The model normalises and converts the fractions, then calculates the mixture properties required by the selected source model (including molecular weight, heat capacities, heat-capacity ratio, and specific gas constant where database inputs are available).

**03. Calculate the source term**

Choose one source-term option:

- **A. Direct mass rate:** supply the gas release rate and release duration.
- **B. Initial gas-orifice estimate:** calculate an initial single-phase ideal-gas discharge rate from upstream conditions and orifice geometry.
- **C. Gas-vessel blowdown:** calculate a time-varying discharge from a rigid, adiabatic, well-mixed ideal-gas vessel.

The source-term calculation provides the gas mass flow as a function of time for the dispersion calculation. The orifice and vessel options do not model liquefied-gas flashing or two-phase discharge.

**04. Calculate the atmospheric wind profile**

- Select the power-law or Monin–Obukhov-type wind-profile option.
- Calculate wind velocity as a function of elevation using the selected stability and surface-roughness inputs.
- Use the resulting wind field in the applicable dispersion calculation.

**05. Select the dispersion engine**

- **A. Gaussian engine:** use the reflected-ground Gaussian plume for a continuous source or finite-duration Gaussian cloud integration for a transient release.
- **B. DEGADISx ground-source engine:** run the external dense-gas model for an eligible, resolved, low-momentum, ground-level area source over flat terrain. The returned downwind profile is mapped onto a 3D display grid for visualisation; this mapping is not a 3D CFD solution.

The Richardson-number calculation is an applicability screen. It does not itself calculate heavy-gas dispersion and does not automatically convert a Gaussian result into a DEGADISx result.

**06. Calculate concentrations**

- Evaluate the concentration field in the downwind (x), crosswind (y), and vertical (z) directions.
- For transient cases, evaluate concentrations at the selected observation times and calculate concentration histories at specified receptors.
- Include ground reflection in the Gaussian formulation.

**07. Perform hazard screening**

- Calculate component concentrations from the specified fixed mixture composition.
- Compare combustible components with available LFL/UFL data and calculate the indicative mixed-fuel limits using Le Chatelier's rule when the required component data are available.
- Compare component concentrations with available AEGL-2, ERPG-2, or IDLH reference values, while respecting each endpoint's meaning and exposure-time basis.
- Treat missing database values as unavailable; do not interpret them as zero or as evidence of no hazard.

**08. Review and export results**

- Review the 3D concentration field, ground/receptor-height contours, vertical cross-section, wind profile, and source-rate or vessel histories as applicable.
- Review receptor concentration histories and peak values.
- Export the available concentration fields, source profiles, and chemical/hazard data to CSV.

**Calculation path in one line:**

`SCENARIO + METEOROLOGY + CHEMICAL MIXTURE -> SOURCE TERM -> WIND PROFILE -> DISPERSION ENGINE -> CONCENTRATION FIELD -> HAZARD SCREENING -> PLOTS / RECEPTORS / CSV`

### 1.2 Dispersion engines

| Engine | Method | Appropriate use | Selection |
|---|---|---|---|
| Gaussian | Reflected-ground Gaussian plume for a steady source; finite-duration Gaussian cloud integration for transient releases. | Passive gas for which gravity does not significantly control cloud motion. | `CONFIG["model_engine"] = "gaussian"` |
| DEGADISx ground source | External DEGADISx package performs a dense-gas downwind-profile calculation. The notebook maps profile variables onto a Cartesian grid for visualisation. | Resolved, low-momentum, ground-level area source over flat, level terrain. | `CONFIG["model_engine"] = "degadisx_ground"` |

The Richardson-number calculation described below is an **applicability screen only**. It does not turn a Gaussian result into a heavy-gas result. The dense-gas engine must be selected explicitly and its source constraints must be satisfied.

### 1.3 Intended applications

- Research, education, and preliminary scenario comparison.
- Examining the effects of release rate, duration, wind speed, stability class, roughness, chemical composition, and receptor position.
- Time-history screening at user-defined receptors.
- First-pass source calculations for direct gas rates, ideal-gas orifice flow, and single-phase ideal-gas vessel blowdown.
- Component-resolved screening against available flammability and toxicological reference values.
- Exploratory dense-gas assessment for a pre-resolved ground-level area source using DEGADISx.

---

## 2. Features and model-selection logic

| Feature group | Available functionality |
|---|---|
| Chemical / mixture | 28 selectable chemical records; single-component or multicomponent mixtures; mole- or mass-fraction inputs; automatic molecular weight, Cp, Cv, gamma, and specific gas constant. |
| Source terms | Direct mass rate; constant-pressure compressible ideal-gas orifice estimate; rigid, adiabatic, well-mixed ideal-gas vessel blowdown. |
| Meteorology | Pasquill stability classes A–F; rural/urban dispersion coefficients; power-law or optional Monin–Obukhov-type wind profile; reference wind and surface roughness. |
| Gaussian dispersion | Steady continuous plume; finite-duration transient cloud; time-varying source segments for vessel blowdown; ground reflection. |
| Computational grid | Nonuniform x-y-z grid; near-source downwind refinement; crosswind centreline clustering; near-ground vertical clustering; exact source and receptor-height insertion. |
| Receptors | User-defined receptors; direct time-dependent concentration histories; selected-time receptor output; direct-versus-grid interpolation comparison. |
| Visualisation | 3D isosurface; ground/receptor-height contours; vertical centreline cross-section; wind profile; release-rate history; vessel pressure and temperature histories. |
| Hazard screening | Component concentrations; indicative Le Chatelier mixed-fuel LFL/UFL; component LFL/UFL reference levels; Fuller diffusivity estimates; available component AEGL-2, ERPG-2 and IDLH endpoints. |
| Dense gas | Optional DEGADISx 0.1.4 integration; steady or transient ground-area source; reconstructed 3D display field; optional receptor histories. |
| Exports | Full 3D snapshot CSV; release-profile CSV; selected-component hazard/database CSV. |

### 2.1 Key selection settings

| Setting | Allowed values | Meaning |
|---|---|---|
| `model_engine` | `gaussian` / `degadisx_ground` | Select passive Gaussian or optional dense-gas engine. |
| `model_mode` | `continuous` / `puff` | Steady plume or finite-duration transient Gaussian cloud. |
| `source_term_mode` | `direct_rate` / `orifice_initial` / `vessel_gas_blowdown` | Select the source-rate model. |
| `mixture.fraction_basis` | `mole` / `mass` | Basis of component fractions supplied by the user. |
| `meteorology.wind_profile_model` | `power_law` / `monin_obukhov` | Select the height-dependent wind profile. |
| `hazard_assessment.toxicity_endpoint` | `AEGL2_60min` / `ERPG2` / `IDLH` | Select the toxicological endpoint family. |

---

## 3. Installation and running in Google Colab

### 3.1 Software and packages

The Gaussian pathway uses Python 3 and the following packages:

- NumPy
- pandas
- SciPy
- Plotly
- IPython display utilities

Google Colab normally provides these packages. The notebook contains an optional installation line if an import is unavailable.

The dense-gas option installs the pinned package `degadisx==0.1.4` when `model_engine` is set to `degadisx_ground`. This requires an internet-enabled Colab runtime. The package publisher describes DEGADISx as a beta Python reimplementation of EPA DEGADIS 2.1 and specifies Python 3.10 or later.

### 3.2 Run the notebook

1. Open [Google Colab](https://colab.research.google.com/).
2. Select **File → Upload notebook** and upload `3D_Gas_Dispersion_Model.ipynb`.
3. Connect the runtime. For DEGADISx, confirm that the runtime has internet access so the pinned package can be installed.
4. Run the notebook from top to bottom using **Runtime → Run all**. Do not skip the chemical database, validation, or function-definition cells.
5. Inspect input-validation messages and warnings, release-rate history, vessel pressure/temperature history (where applicable), wind profile, grid diagnostics, concentration plots, receptor tables, and hazard outputs.
6. For a new scenario, edit **Section 2 — INPUT CONFIGURATION** and run the notebook from the top again.

> **Reproducibility:** Save a separate copy of the notebook for each scenario. Record the selected chemicals and fractions, source settings, meteorology, engine, observation time, manual overrides, notebook revision/date, and exported results.

---

## 4. Quick-start configuration

Normal scenario edits should be made in **Section 2 — INPUT CONFIGURATION**. The following example uses a 70 mol% propane / 30 mol% n-butane vapour mixture with a direct finite-duration gas release.

```python
CONFIG["model_engine"] = "gaussian"
CONFIG["model_mode"] = "puff"

CONFIG["source"].update({
    "source_term_mode": "direct_rate",
    "mass_rate_kg_s": 0.05,
    "release_duration_s": 120.0,
    "release_height_m": 1.0,
    "source_y_m": 0.0,
})

CONFIG["mixture"] = {
    "fraction_basis": "mole",
    "components": [
        {"chemical": "Propane (C3H8)", "fraction": 0.70},
        {"chemical": "n-Butane (n-C4H10)", "fraction": 0.30},
    ],
}

CONFIG["meteorology"].update({
    "wind_speed_reference_m_s": 3.0,
    "wind_reference_height_m": 10.0,
    "stability_class": "D",
    "terrain": "Rural",
    "wind_profile_model": "power_law",
})

CONFIG["puff"].update({
    "observation_time_s": 60.0,
    "history_end_s": 600.0,
    "history_step_s": 5.0,
})
```

### 4.1 Select chemicals and composition

Choose exact chemical names from the embedded database. The names are case-sensitive. Each species should appear once, and the supplied fractions must be nonnegative and sum to **1.0** or **100**. Components with zero fraction are ignored.

For mass fractions, use:

```python
CONFIG["mixture"]["fraction_basis"] = "mass"
```

The notebook converts the supplied fraction basis into both mole and mass fractions before calculating mixture properties.

### 4.2 Choose the source model

- `direct_rate`: user supplies gas mass rate in kg/s.
- `orifice_initial`: calculates an initial single-phase ideal-gas orifice rate at constant upstream pressure for the configured release duration.
- `vessel_gas_blowdown`: calculates a time-dependent release from a finite gas inventory in a rigid, adiabatic, ideal-gas vessel. Use `model_mode = "puff"` for vessel blowdown.

### 4.3 Choose the dispersion engine

Use `gaussian` for passive-gas plume/puff calculations. Use `degadisx_ground` only if the source already represents a resolved, low-momentum, ground-level area source suitable for the DEGADISx ground-source model. The dense-gas route has stricter constraints described in [Section 5.7](#57-dense-gas-settings-and-degadisx-source-constraints).

---

## 5. Input configuration reference

All values below describe the notebook's embedded defaults. Units are SI unless stated otherwise. Confirm the actual values in Section 2 of the notebook before running a scenario.

### 5.1 Source settings


| Input key | Default | Unit / meaning |
| --- | --- | --- |
| source.source_term_mode | direct_rate | Direct input, constant-pressure orifice, or ideal-gas vessel blowdown. |
| source.mass_rate_kg_s | 0.05 | kg/s; release rate used by direct_rate, or initial rate for orifice_initial. |
| source.release_duration_s | 120 | s; duration for direct_rate and orifice_initial in transient mode. |
| source.release_height_m | 1.0 | m; source height. Dense-gas ground-area mode requires 0 m. |
| source.source_y_m | 0.0 | m; crosswind offset of source in local plume coordinates. |
| source.orifice.upstream_pressure_bar_abs | 5.0 | bar absolute; upstream pressure for constant-pressure gas-orifice estimate. |
| source.orifice.hole_diameter_mm | 5.0 | mm; circular orifice diameter. This is not the DEGADISx source radius. |
| source.orifice.discharge_coefficient | 0.80 | Dimensionless Cd. |
| source.orifice.gas_temperature_K | 298.15 | K; upstream gas temperature. |
| source.vessel.vessel_volume_m3 | 0.10 | m³; rigid vessel volume. |
| source.vessel.initial_pressure_bar_abs | 8.0 | bar absolute; vessel initial pressure. |
| source.vessel.initial_temperature_K | 298.15 | K; initial bulk temperature. |
| source.vessel.hole_diameter_mm | 3.0 | mm; circular leak opening. |
| source.vessel.discharge_coefficient | 0.80 | Dimensionless Cd for the vessel leak orifice. |
| source.vessel.heat_capacity_ratio_gamma_override | None | Optional vessel-flow gamma override; None uses calculated mixture gamma. |
| source.vessel.maximum_release_duration_s | 1800 | s; upper integration duration if pressure has not approached ambient. |
| source.vessel.profile_output_step_s | 1.0 | s; output spacing for the blowdown profile. |
| source.orifice.heat_capacity_ratio_gamma_override | None | Optional constant-pressure orifice gamma override; None uses calculated mixture gamma. |


### 5.2 Meteorology


| Input key | Default | Meaning / guidance |
| --- | --- | --- |
| wind_speed_reference_m_s | 3.0 | m/s; measured wind speed at the selected reference height. |
| wind_reference_height_m | 10.0 | m; height at which the reference wind speed applies. Must exceed roughness length. |
| wind_direction_from_deg | 0.0 | Degrees meteorological “from”: 0° N, 90° E. Used in geographic offsets in the full-field CSV; model plots are wind-aligned. |
| temperature_K | 298.15 | K; ambient temperature for ideal-gas density and concentration-unit conversion. |
| pressure_Pa | 101325 | Pa absolute; ambient pressure. |
| relative_humidity_pct | 50 | %; passed to DEGADISx deck. |
| air_molecular_weight_g_mol | 28.965 | g/mol; ambient-air molecular weight used in ideal-gas density and Richardson screening. |
| stability_class | D | Pasquill class A-F. A-C unstable, D neutral, E-F stable. |
| terrain | Rural | Rural or Urban; selects the vertical Briggs dispersion coefficient set. |
| roughness_length_m | 0.03 | m; surface roughness length used in the wind-profile/stability parameterisation. |
| wind_profile_model | power_law | power_law or monin_obukhov. |
| power_law_n_override | None | Optional site-calibrated exponent. If None, class-specific defaults are used. |


### 5.3 Computational domain and receptors


| Setting | Default | Meaning |
| --- | --- | --- |
| domain.x_max_m | 1200 | m; downwind computational extent. |
| domain.y_half_width_m | 200 | m; crosswind domain extends from -value to +value. |
| domain.z_max_m | 150 | m; maximum height. |
| domain.nx / ny / nz | 180 / 91 / 65 | Base counts for nonuniform x, y and z grids; exact source/receptor coordinates may be inserted. |
| domain.x_refinement_end_m | 30 | m; near-source linear-grid extent before geometric stretching. |
| domain.y_stretch_factor | 2.5 | Controls clustering around the y centreline. |
| domain.z_stretch_power | 1.6 | Controls clustering near ground level. |
| domain.receptor_height_m | 1.5 | m; inserted into z grid for ground/receptor-height contours. |
| RECEPTORS | 3 examples | List of receptor names and x_m, y_m, z_m in local wind-aligned coordinates. Edit or add receptors. |


The default receptor list is:

| Receptor | x (m) | y (m) | z (m) |
|---|---:|---:|---:|
| Receptor 1 | 100 | 0 | 1.5 |
| Receptor 2 | 250 | 25 | 1.5 |
| Receptor 3 | 500 | −25 | 1.5 |

Receptors use local wind-aligned coordinates. Add or edit entries in the `RECEPTORS` list.

### 5.4 Transient and visualisation settings


| Setting | Default | Meaning |
| --- | --- | --- |
| puff.release_interval_s | 1.0 | s; time-bin used to approximate a variable blowdown-rate history. A constant direct/orifice release is integrated as one exact source interval. |
| puff.observation_time_s | 60 | s; snapshot time for Gaussian transient field. Default is during the 120 s release. |
| puff.history_end_s | 600 | s; end of receptor history. Increase to at least 3600 s for AEGL-2/ERPG-2 moving-window screening or 1800 s for IDLH window screening. |
| puff.history_step_s | 5 | s; time step for receptor history. |
| plot.isosurface_level_ppm | 10 | ppmv; 3D display threshold only, not a safety criterion. |
| plot.ground_contour_level_ppm | 10 | ppmv; contour reference value. |
| plot.contour_scale | log | log or linear concentration colour scale. |
| plot.log_floor_ppm | 0.01 | ppmv; lower colour display limit only. It does not change the physical concentration field. |
| plot.show_3d_isosurface | True | Turn off for very large fields if browser rendering is slow. |


The isosurface and contour levels are display settings. They are not safety criteria unless intentionally set equal to a relevant reference endpoint. For logarithmic plots, `log_floor_ppm` controls the minimum display value only; it does not change calculated concentrations.

### 5.5 Hazard-assessment settings


| Setting | Default | Meaning / interpretation |
| --- | --- | --- |
| hazard_assessment.plot_hazard_contours | True | Enable/disable component-specific hazard contour calculations and plots. |
| hazard_assessment.toxicity_endpoint | AEGL2_60min | Select AEGL-2 (60 min), ERPG-2 or IDLH where a numerical component endpoint exists. |
| hazard_assessment.flammability_rule | LeChatelier | Selects the approximate Le Chatelier fuel-mixture LFL/UFL rule. Missing required component limits disable the blend calculation. |
| hazard_assessment.snapshot_contours_are_instantaneous | True | The plotted concentration snapshot is compared instantaneously with the selected reference; it is not a duration-averaged exposure map. |


### 5.6 Dense-gas settings


| Setting | Default | Meaning / important restriction |
| --- | --- | --- |
| heavy_gas.source_radius_m | 1.0 m | Equivalent radius of the resolved ground-level gas source plane. Not an orifice/hole diameter. |
| heavy_gas.source_temperature_K | 298.15 K | Resolved source temperature. |
| heavy_gas.averaging_time_s | 600 s | Averaging-time input to DEGADISx. A 10 s value may be used for flammable-cloud averaging only where justified. |
| heavy_gas.mole_fraction_cutoff | 1e-6 | Numerical contaminant concentration cutoff. |
| heavy_gas.snapshot_selection | peak_centerline | peak_centerline or closest_to_time. |
| heavy_gas.snapshot_time_s | 120 s | Requested time when closest_to_time is selected. |
| heavy_gas.nx / ny / nz | 100 / 81 / 51 | Resolution of the reconstructed visualisation grid, not the native DEGADISx solver grid. |
| heavy_gas.y_half_width_m | 200 m | Half-width of the post-processing crosswind field. |
| heavy_gas.z_max_m | 80 m | Height of the post-processing field. |
| heavy_gas.plot_isosurface_level_ppm | 100 ppmv | Display threshold, not a flammability or toxic endpoint. |
| heavy_gas.run_receptor_histories | True | Run DEGADISx transient receptor histories when transient mode is selected. |


> **DEGADISx source requirements:** `model_engine = "degadisx_ground"` requires `source_term_mode = "direct_rate"`, `release_height_m = 0`, `source_y_m = 0`, a positive `source_radius_m`, and a positive direct mass rate. `source_radius_m` is the equivalent radius of the resolved ground-level source plane. It is **not** the leak-hole/orifice diameter.
>
> Resolve the high-pressure jet, flashing, droplet transport, aerosol rainout, and any pool evaporation upstream. Enter only the resulting resolved ground-level vapour-source rate and area/radius. This notebook does not perform that bridge. The dense-gas route is intended for flat, level terrain and low-momentum area sources.

### 5.7 Export settings


| Setting | Default | What is written |
| --- | --- | --- |
| export.write_full_3d_csv | False | Full concentration snapshot at all 3D grid points, including local x-y-z and east/north offsets, mixture kg/m³, mg/m³ and ppmv, plus component mass concentration and ppmv columns. Filename key: csv_filename. |
| export.write_release_profile_csv | False | Release time and mass-rate profile; vessel blowdown also includes tank inventory, pressure, temperature and flow regime. Filename key: release_profile_filename. |
| export.write_hazard_database_csv | False | Selected-component hazard and transport-property table with data-source links. Filename key: hazard_database_filename. |


---

## 6. Model equations and methods

This section documents the equations implemented in the notebook. Equations are written in SI units unless stated otherwise. “ALOHA-based” means a parameterisation or formulation is taken from or inspired by the supplied ALOHA 5.4.4 technical documentation; it does **not** mean that the notebook duplicates ALOHA source code or all of its numerical algorithms.

### 6.1 Ideal-gas and gas-mixture properties

The universal gas constant is:

```math
R_u = 8.314462618\;\mathrm{J\,mol^{-1}\,K^{-1}}
```

**Ideal-gas density**

```math
\rho = \frac{PM}{R_uT}
```

where `ρ` is density (kg/m³), `P` is absolute pressure (Pa), `M` is molecular weight (kg/mol), and `T` is absolute temperature (K). The same equation is applied to ambient air using the configured air molecular weight. Ideal-gas behaviour is assumed.

**Specific gas constant**

```math
R_s = \frac{R_u}{M_{mix}}
```

where `M_mix` is the mixture molecular weight in kg/mol and `R_s` is in J/(kg K).

**Mole-fraction input**

```math
M_{mix}=\sum_i y_iM_i
```

```math
w_i=\frac{y_iM_i}{M_{mix}}
```

Here `y_i` and `w_i` are the mole and mass fractions of component `i`.

**Mass-fraction input**

```math
M_{mix}=\left[\sum_i\frac{w_i}{M_i}\right]^{-1}
```

```math
y_i=\frac{w_i/M_i}{\sum_j w_j/M_j}
```

Fractions are normalised to a sum of unity. The mixture is assumed ideal and composition remains fixed during release and dispersion; reaction, vapour-liquid fractionation, and differential species loss are not modelled.

**Mixture heat capacity and heat-capacity ratio**

```math
C_{p,mix}=\sum_i y_iC_{p,i}
```

```math
C_{v,mix}=C_{p,mix}-R_u,\qquad \gamma=\frac{C_{p,mix}}{C_{v,mix}}
```

The database Cp values are gas-phase values at or near 298.15 K and are treated as constants. The calculated `γ` is used in ideal-gas orifice flow and vessel blowdown unless an override is supplied. This approximation becomes less reliable when the gas temperature changes substantially or non-ideal effects matter.

### 6.2 Atmospheric wind profiles

Wind speed changes with height because of surface drag and atmospheric stability. The entered wind speed is recovered at the specified reference height.

#### 6.2.1 Power-law profile

```math
U(z)=U_{ref}\left(\frac{z}{z_{ref}}\right)^n
```

`U(z)` is wind speed at height `z`, `U_ref` is the measured speed at reference height `z_ref`, and `n` is a dimensionless stability exponent. The embedded defaults are typical ALOHA-documented values for roughness length `z_0=0.15` m; use site-calibrated values where the actual surface differs materially.

| Pasquill class | Description | Default exponent `n` |
|---|---|---:|
| A | Extremely unstable | 0.109 |
| B | Moderately unstable | 0.112 |
| C | Slightly unstable | 0.120 |
| D | Neutral | 0.142 |
| E | Slightly stable | 0.203 |
| F | Stable | 0.253 |

#### 6.2.2 Monin–Obukhov-type profile

The alternative profile is based on a logarithmic surface-layer wind law with a stability correction. The Obukhov length is estimated from roughness and Pasquill class:

| Class | Obukhov length `L` (m) |
|---|---|
| A | `-11.4z_0^0.10` |
| B | `-26.0z_0^0.17` |
| C | `-123z_0^0.30` |
| D | `∞` (neutral) |
| E | `+123z_0^0.30` |
| F | `+26.0z_0^0.17` |

Define the dimensionless height:

```math
\zeta=\frac{z}{L}
```

The implemented momentum stability function is:

```math
\psi_m(\zeta)=0,\quad \zeta=0
```

```math
\psi_m(\zeta)=2\ln\left(\frac{1+a}{2}\right)+\ln\left(\frac{1+a^2}{2}\right)-2\tan^{-1}(a)+\frac{\pi}{2},\quad \zeta<0
```

```math
a=(1-16\zeta)^{1/4}
```

```math
\psi_m(\zeta)=-5\zeta,\quad \zeta>0
```

The code uses the Businger–Dyer-style sign convention in which the wind profile contains `-ψ_m`. This is a related surface-layer parameterisation, not a bit-for-bit reproduction of ALOHA's exact equations.

The friction velocity is inferred from the reference wind:

```math
u_* = \frac{\kappa U_{ref}}{\ln[(z_{ref}+z_0)/z_0]-\psi_m(z_{ref}/L)}
```

The resulting wind profile is:

```math
U(z)=\frac{u_*}{\kappa}\max\left\{\ln\left(\frac{z+z_0}{z_0}\right)-\psi_m(z/L),\;0\right\}
```

where `κ=0.40` is the von Kármán constant. A separate neutral-log friction-velocity estimate is used by the Richardson-number screen.

#### 6.2.3 Effective plume wind speed

The Gaussian solver uses a vertical concentration-weighted wind speed at each downwind station rather than wind speed at one arbitrary elevation. Define:

```math
W(z;x)=\exp\left[-\frac{(z-H)^2}{2\sigma_z^2}\right]+\exp\left[-\frac{(z+H)^2}{2\sigma_z^2}\right]
```

```math
U_{eff}(x)=\frac{\int_0^{z_{cap}}U(z)W(z;x)\,dz}{\int_0^{z_{cap}}W(z;x)\,dz}
```

`H` is release height and `σ_z` is vertical dispersion width. The code evaluates the integrals using trapezoidal quadrature over 360 vertical samples, with `z_cap=max(250 m, H+8max(σ_z))`. A lower transport-speed limit of 0.05 m/s avoids division by zero in the numerical calculations.

### 6.3 Briggs / ALOHA dispersion parameters

The Gaussian plume widths depend on downwind distance, stability class and roughness. The coefficients below are transcribed from Table 13 in the supplied ALOHA 5.4.4 documentation, which presents Briggs-based empirical relationships.

```math
\sigma_x(x)=s_{x1}x^{s_{x2}}
```

```math
\sigma_y(x)=\frac{s_{y1}x}{\sqrt{1+s_{y2}x}}
```

```math
\sigma_z(x)=s_{z1}x(1+s_{z2}x)^{s_{z3}}
```

All `σ` values are standard deviations in metres. The code evaluates the point-source equations at a minimum downwind distance of 1 m to regularise the source-plane calculation. This is a numerical treatment, not a physical near-field or nozzle-jet model.

**Coefficients for `σ_x` and `σ_y`**


| Stability | sₓ₁ | sₓ₂ | sᵧ₁ | sᵧ₂ |
| --- | --- | --- | --- | --- |
| A | 0.02 | 1.22 | 0.22 | 0.0001 |
| B | 0.02 | 1.22 | 0.16 | 0.0001 |
| C | 0.02 | 1.22 | 0.11 | 0.0001 |
| D | 0.04 | 1.14 | 0.08 | 0.0001 |
| E | 0.17 | 0.97 | 0.06 | 0.0001 |
| F | 0.17 | 0.97 | 0.04 | 0.0001 |


**Vertical dispersion coefficients — rural terrain**


| Stability | s_z₁ | s_z₂ | s_z₃ |
| --- | --- | --- | --- |
| A | 0.2 | 0 | 0 |
| B | 0.12 | 0 | 0 |
| C | 0.08 | 0.0002 | -0.5 |
| D | 0.06 | 0.0015 | -0.5 |
| E | 0.03 | 0.0003 | -1 |
| F | 0.016 | 0.0003 | -1 |


**Vertical dispersion coefficients — urban terrain**


| Stability | s_z₁ | s_z₂ | s_z₃ |
| --- | --- | --- | --- |
| A | 0.24 | 0.001 | 0.5 |
| B | 0.24 | 0.001 | 0.5 |
| C | 0.2 | 0 | 0 |
| D | 0.14 | 0.0003 | -0.5 |
| E | 0.08 | 0.0015 | -0.5 |
| F | 0.08 | 0.0015 | -0.5 |


The implementation uses the rural `s_y` coefficients for both terrain classes, following the ALOHA discussion of differing averaging times for the source datasets. The Rural/Urban selection changes the `σ_z` coefficient set.

### 6.4 Continuous Gaussian plume

For a steady point source with mass release rate `Q`, the code uses a reflected-ground Gaussian plume:

```math
\begin{aligned}
C(x,y,z)=&\;\frac{Q}{2\pi U_{eff}\sigma_y\sigma_z}
\exp\left[-\frac{(y-y_s)^2}{2\sigma_y^2}\right]\\
&\times\left\{\exp\left[-\frac{(z-H)^2}{2\sigma_z^2}\right]+\exp\left[-\frac{(z+H)^2}{2\sigma_z^2}\right]\right\}
\end{aligned}
```

`C` is total mixture mass concentration (kg/m³), `Q` is mass release rate (kg/s), `U_eff` is m/s, and `y_s` and `H` are the source crosswind offset and release height. The second vertical exponential is an image source below ground; it represents ground reflection in the idealised solution. The equation is intended for passive pollutants and does not capture dense-gas slumping or source-jet momentum.

### 6.5 Finite-duration transient Gaussian cloud

The transient calculation integrates the longitudinal Gaussian distribution over each finite source interval. It is inspired by the short-duration Gaussian-cloud formulation in Section 4.3 of the ALOHA technical documentation and by Palazzi et al. (1982). It is not a direct reimplementation of ALOHA's proprietary code.

For source segment `i` with average mass rate `\bar q_i` released from `t_0,i` to `t_1,i` and observed at time `t`:

```math
t_{1,eff}=\min(t_{1,i},t)
```

```math
a_{latest}=\max(t-t_{1,eff},0),\qquad a_{earliest}=\max(t-t_{0,i},0)
```

```math
G_{x,i}=\frac{1}{2}\left\{\operatorname{erf}\left[\frac{x-U_{eff}a_{latest}}{\sqrt{2}\sigma_x}\right]-\operatorname{erf}\left[\frac{x-U_{eff}a_{earliest}}{\sqrt{2}\sigma_x}\right]\right\}
```

The crosswind and vertical kernels are:

```math
g_y=\frac{1}{\sqrt{2\pi}\sigma_y}\exp\left[-\frac{(y-y_s)^2}{2\sigma_y^2}\right]
```

```math
g_z=\frac{1}{\sqrt{2\pi}\sigma_z}\left\{\exp\left[-\frac{(z-H)^2}{2\sigma_z^2}\right]+\exp\left[-\frac{(z+H)^2}{2\sigma_z^2}\right]\right\}
```

The total transient field is the sum of source-segment contributions:

```math
C(x,y,z,t)=\sum_i\left[\frac{\bar q_i}{U_{eff}}G_{x,i}g_yg_z\right]
```

`erf` is the Gauss error function. The age terms describe the youngest and oldest gas emitted during a segment. If observation occurs during active release, the segment is clipped at the observation time and the source plane can have non-zero concentration. After shutoff, the cloud is expected to advect away from the source.

For `direct_rate` and `orifice_initial`, the transient source is an exact finite-duration top-hat interval. For vessel blowdown, the time-varying mass-rate history is divided into intervals. The released mass in an interval is approximated by trapezoidal integration:

```math
\Delta m_i\approx\frac{1}{2}(\dot m_i+\dot m_{i+1})\Delta t_i,\qquad \bar q_i=\frac{\Delta m_i}{\Delta t_i}
```

Smaller release intervals resolve changes in blowdown rate more closely but increase the number of source segments.

### 6.6 Concentration units and component concentrations

The model converts mixture mass concentration to ppmv using the ideal-gas relation:

```math
C_{ppmv}=\left(\frac{C_{kg/m^3}}{M_{kg/mol}}\right)\left(\frac{R_uT}{P}\right)10^6
```

For a fixed-composition mixture, component fields are estimated as:

```math
C_{i,ppmv}(x,y,z,t)=y_iC_{mix,ppmv}(x,y,z,t)
```

```math
C_{i,kg/m^3}(x,y,z,t)=w_iC_{mix,kg/m^3}(x,y,z,t)
```

This partition is valid only under the assumption that composition remains constant everywhere and at all times. Differential component diffusion, condensation, evaporation, reaction, and preferential sorption are not represented.

### 6.7 Gas-release equations

#### 6.7.1 Direct mass rate

```math
\dot m(t)=\dot m_{user},\qquad 0\le t\le t_{release}
```

The mass rate is supplied directly in kg/s. In transient mode it is a top-hat source over `release_duration_s`; in continuous mode it supplies the steady plume.

#### 6.7.2 Compressible ideal-gas orifice

For a single-phase perfect gas through a circular orifice:

```math
A=\frac{\pi d^2}{4},\qquad R_s=\frac{R_u}{M_{mix}}
```

```math
r=\frac{P_{down}}{P_{up}},\qquad r_{crit}=\left[\frac{2}{\gamma+1}\right]^{\gamma/(\gamma-1)}
```

If `r≤ r_crit`, flow is choked:

```math
\dot m=C_dAP_{up}\sqrt{\frac{\gamma}{R_sT_{up}}}\left[\frac{2}{\gamma+1}\right]^{\frac{\gamma+1}{2(\gamma-1)}}
```

If `r>r_crit`, flow is unchoked:

```math
\dot m=C_dAP_{up}\sqrt{\frac{2\gamma}{R_sT_{up}(\gamma-1)}\left[r^{2/\gamma}-r^{(\gamma+1)/\gamma}\right]}
```

Here `d` is the hole diameter, `C_d` is the discharge coefficient, `γ` is the heat-capacity ratio, `T_up` is upstream temperature, and pressures are absolute. `orifice_initial` holds the initial upstream pressure constant for the configured release duration. It is not a tank blowdown calculation.

#### 6.7.3 Adiabatic rigid-vessel blowdown

Initial vessel mass is estimated using:

```math
m_0=\frac{P_0V}{R_sT_0}
```

For a rigid, adiabatic, well-mixed ideal-gas vessel with constant `γ` and fixed composition:

```math
\frac{P}{P_0}=\left(\frac{m}{m_0}\right)^\gamma
```

```math
\frac{T}{T_0}=\left(\frac{m}{m_0}\right)^{\gamma-1}
```

The compressible-orifice equation is evaluated at the current vessel pressure and temperature. The vessel mass balance is:

```math
\frac{dm}{dt}=-\dot m_{orifice}\big(P(m),P_{ambient},T(m)\big)
```

The ODE is integrated with SciPy `solve_ivp` (relative tolerance `10^-7`). Integration stops when vessel pressure approaches 1.0001 times ambient pressure or when the configured maximum blowdown time is reached. The release profile includes time, tank mass, absolute pressure, temperature, discharge rate, and choked/unchoked regime.

### 6.8 Richardson-number applicability screen

The model calculates an ALOHA-style critical Richardson number as an indicator of whether gravitational effects may be important:

```math
g'=g\frac{\rho_g-\rho_{air}}{\rho_{air}}
```

```math
H^*=\frac{\dot m\pi}{4\rho_gU_{10}}
```

```math
Ri_c^*=\frac{H^*g'}{u_*^2}
```

Here `ρ_g` is ideal-gas release-mixture density at ambient conditions, `ρ_air` is ambient-air density, `g=9.80665` m/s², `U_10` is wind speed at 10 m from the selected profile, and `u_*` is the neutral-log friction-velocity estimate used in this screen. The notebook indicates passive-Gaussian behaviour when `Ri_c^*<1` and warns that dense-gas effects may be important when `Ri_c^*≥1`. This simple screen does not include source-temperature or real-gas density effects and is **not a dense-gas solver**.

### 6.9 Fuller–Schettler–Giddings binary diffusivity

Where a Fuller diffusion volume is available, binary molecular diffusivity of component `i` in air is estimated using:

```math
M_{ij}=\frac{2}{1/M_i+1/M_j}
```

```math
D_{ij}=\frac{0.00143T^{1.75}}{P\sqrt{M_{ij}}\left[V_i^{1/3}+V_j^{1/3}\right]^2}
```

Implemented units: `T` in K, `P` in atm, `M` in g/mol, and Fuller volumes `V` in cm³/mol; the result is cm²/s and is multiplied by `10^-4` to report m²/s. Air is represented by MW = 28.965 g/mol and Fuller volume = 19.7 cm³/mol. Some volumes are tabulated, while others are group/atomic-increment estimates. They are not all measured pair-specific binary diffusivities.

Molecular diffusivity is reported as a component property. It is **not added to** the empirical Briggs `σ_y` or `σ_z` values, because those atmospheric dispersion parameters already represent turbulent spreading; adding molecular diffusivity again could double-count mixing.

### 6.10 Flammability screening

For combustible components with numerical LFL and UFL values, the notebook applies the Le Chatelier mixing rule. First normalise the fuel-component mole fractions:

```math
y_{i,fuel}=\frac{y_i}{\sum_{j\in fuel}y_j}
```

Then calculate the fuel-submixture limits:

```math
LFL_{mix}=\left[\sum_{i\in fuel}\frac{y_{i,fuel}}{LFL_i}\right]^{-1}
```

```math
UFL_{mix}=\left[\sum_{i\in fuel}\frac{y_{i,fuel}}{UFL_i}\right]^{-1}
```

The corresponding threshold on total released-mixture concentration is estimated by:

```math
LFL_{total,cloud}=\frac{LFL_{mix}}{y_{fuel}},\qquad UFL_{total,cloud}=\frac{UFL_{mix}}{y_{fuel}}
```

LFL/UFL values are stored as vol% in air; multiplication by 10,000 converts vol% to ppmv. `y_fuel` is the total mole fraction of combustible components in the released mixture. The rule is an approximate screening method: it does not represent oxygen enrichment, inerting, temperature/pressure corrections, ignition energy, chemical interactions, or non-ideal flammability behaviour.

If a component is identified as flammable/pyrophoric but lacks a numeric LFL or UFL, the blend calculation is disabled rather than silently omitting that species. For example, phosphine has a listed LFL but no numeric UFL in the curated source. Oxygen is treated as an oxidiser, not a fuel component; oxygen-enrichment effects are not calculated.

### 6.11 Toxicity endpoints and receptor exposure screening

The notebook provides component-specific reference endpoints where database values are available:

- AEGL-2: 60-minute exposure endpoint;
- ERPG-2: 1-hour emergency-planning guideline;
- IDLH: 30-minute escape-planning reference.

Component concentration is estimated as:

```math
C_i(t)=y_iC_{mix}(t)
```

For receptor histories, a moving average is calculated over the selected exposure window:

```math
\overline C_{i,\Delta t}(t_k)=\frac{1}{N}\sum_{j=0}^{N-1}C_i(t_{k-j}),\qquad N=\frac{\Delta t}{\Delta t_{sample}}+1
```

The AEGL-2 and ERPG-2 windows are 3,600 s; the IDLH screening window is 1,800 s. A result is only assessed when the simulated receptor history covers the full window with sufficient samples. Otherwise, the status is `NOT ASSESSED`. Snapshot hazard contours are instantaneous spatial comparisons, not duration-averaged toxicity maps. The model does not apply an additive mixture-toxicity formula.

IDLH records whose source basis is 10% LEL (including methane, ethane, propane, n-butane and 1,3-butadiene in the embedded data) are excluded from toxic-endpoint contour screening. Those values are not interpreted as acute toxicological thresholds.

### 6.12 Nonuniform grid and display scaling

The grid is stretched to allocate more points near the source, plume centreline, and ground. Base coordinate transforms are:

```math
x_{near}=\operatorname{linspace}(0,x_{ref},n_{near}),\qquad x_{far}=\operatorname{geomspace}(x_{ref},x_{max},n_{far})
```

```math
y(\eta)=y_{half}\frac{\tanh(\beta\eta)}{\tanh(\beta)},\qquad \eta\in[-1,1]
```

```math
z(\eta_z)=z_{max}\eta_z^p,\qquad \eta_z\in[0,1]
```

The exact x = 0 source plane, source crosswind coordinate, release height, ground height and receptor height are inserted into the coordinate arrays. A grid diagnostic compares local spacing with dispersion widths near 10 m downwind. This diagnostic is not proof of grid convergence; perform explicit coarse/medium/fine comparisons before relying on receptor or contour values.

For full-field export, local wind-aligned coordinates are rotated to east/north offsets. With `θ_to=(θ_from+180^°)mod360^°`:

```math
East=x\sin\theta_{to}+y\cos\theta_{to}
```

```math
North=x\cos\theta_{to}-y\sin\theta_{to}
```

This is an export-coordinate transformation only. It does not modify dispersion. Logarithmic display applies `log_10(max(C_ppmv,C_floor))` to the plotted values only; it does not change the physical concentrations.

### 6.13 DEGADISx dense-gas engine and reconstructed 3D display

When `model_engine = "degadisx_ground"`, the notebook writes a DEGADISx input deck and calls the external package's `run_steady` or `run_transient` function (package version 0.1.4). The DEGADIS calculations are performed by the external DEGADISx package. The original DEGADIS model is intended for ground-level dense gas/aerosol area sources with negligible release momentum over flat, level terrain; it represents gravity-driven flow and entrainment.

The package reports downwind profile variables including centreline contaminant mass concentration `c_c(x)`, core half-width `b(x)`, lateral scale `s_y(x)`, vertical scale `s_z(x)`, and vertical profile exponent `α`. The notebook interpolates these outputs onto x-grid stations and reconstructs a 3D display field:

```math
\alpha_1=\alpha+1
```

```math
f_y(x,y)=1,\qquad |y|\le b(x)
```

```math
f_y(x,y)=\exp\left[-\left(\frac{|y|-b(x)}{s_y(x)}\right)^2\right],\qquad |y|>b(x)
```

```math
f_z(x,z)=\exp\left[-\left(\frac{z}{s_z(x)}\right)^{\alpha_1}\right]
```

```math
C_{c,3D}(x,y,z)=c_c(x)f_y(x,y)f_z(x,z)
```

This is a post-processing approximation for displaying and interrogating the downwind profile. It is **not** a native 3D DEGADIS field, a CFD calculation, or a separate three-dimensional dense-gas solver.

The profile-derived contaminant mass concentration is converted to contaminant mole fraction using the DEGADISx thermodynamic table. The notebook interpolates mixture density from the package table:

```math
\rho=\operatorname{interp}(c_c;\;table.c_c,table.\rho)
```

```math
w_c=\frac{c_c}{\rho}
```

```math
w_a=\frac{1-(1+h_s)w_c}{1+h}
```

```math
w_w=1-w_a-w_c
```

```math
M_{mix}=\left[\frac{w_c}{MW_c}+\frac{w_a}{28.9647}+\frac{w_w}{18.01528}\right]^{-1}
```

```math
y_c=\left(\frac{M_{mix}}{MW_c}\right)w_c,\qquad C_{c,ppmv}=10^6y_c
```

Here `c_c` is contaminant mass concentration (kg/m³), `ρ` is total mixture density from the package table, `w_c`, `w_a`, and `w_w` are contaminant, dry-air, and water mass fractions, `h_s` and `h` are source and ambient humidity variables, and `MW_c` is released-gas molecular weight. The 28.9647 and 18.01528 g/mol terms are the air and water molecular weights used by the package conversion. This conversion is distinct from the Gaussian-pathway ideal-gas ppmv conversion.

---

## 7. Outputs and interpretation

| Output | What it shows | Interpretation |
|---|---|---|
| Input validation and mixture summary | Selected components, fractions, mixture MW, Cp, gamma, specific gas constant, endpoint availability. | Confirm composition and data-source notes before interpreting results. |
| Wind profile | Wind speed versus height and reference-height marker. | The profile should recover the entered wind speed at the reference height. Compare with site measurements. |
| Release history | Mass rate versus time; vessel pressure and temperature for blowdown. | Check pressure decrease, flow regime, release termination and whether the selected source model is appropriate. |
| 3D concentration field | Isosurface at selected ppmv display threshold. | Isosurface is a selected concentration level, not necessarily a hazard boundary. |
| Ground/receptor-height contour | Horizontal concentration map at the configured receptor height. | Log scale helps show dilute plume periphery; inspect actual hover values and colour floor. |
| Vertical centreline cross-section | Downwind distance versus elevation. | Shows vertical plume distribution and ground interaction under the model assumptions. |
| Receptor histories | Direct concentration versus time at each receptor. | Use arrival time, peak concentration, and decay. Direct receptor evaluation avoids coarse-grid interpolation. |
| Receptor snapshot table | Direct calculation and grid-interpolated concentration. | Large differences indicate insufficient spatial resolution or receptor coordinates outside the grid. |
| Hazard overlays | Mixed-fuel LFL/UFL and component toxicity endpoint references. | Missing values are not plotted. Toxicity contours are component-specific and instantaneous at the selected snapshot. |
| DEGADISx reconstruction | 3D display from downwind profile quantities. | Treat as reconstructed display only; the underlying calculations are performed by DEGADISx. |

### 7.1 Source-connected plume interpretation

The computational grid includes x = 0, source crosswind coordinate, and release height. During an active finite-duration release, the transient integration can produce concentration at the source plane, subject to the selected isosurface threshold and the 1 m sigma regularisation. After shutoff, the cloud is expected to move downwind and leave a source-to-cloud gap. A visible gap does not automatically prove a code error: inspect source-plane concentration and the display threshold.

### 7.2 CSV exports

- **Full 3D snapshot CSV:** local x/y/z; east/north offsets; mixture concentration in kg/m³, mg/m³ and ppmv; component-specific mass and mole concentration columns. This is the full concentration snapshot grid, not a continuous CFD data file.
- **Release profile CSV:** time and mass rate; vessel inventory, pressure, temperature and choked/unchoked regime for vessel blowdown.
- **Selected hazard database CSV:** selected component fractions, LFL/UFL, diffusivity, toxicity endpoints, notes, and data-source links.

---

## 8. Chemical database

The embedded library contains **28 gas-phase component records**. Thermophysical data include formula, CAS number, molecular weight and ideal-gas heat capacity near 298.15 K. The hazard/transport library includes LFL/UFL where sourced, Fuller diffusion volume, estimated gas-in-air diffusivity, AEGL-2 (60 min), ERPG-2 and IDLH where available.

A blank (`—`) means that no numeric value is populated in this curated database. It does **not** mean zero, nonhazardous, or safe. Chemical-specific source URLs are embedded in the notebook database dictionaries; consult those records for the exact underlying page.

### 8.1 Thermophysical and diffusivity database

The diffusivity column below is $`D_{i-air}\times10^5`$ in m²/s, calculated using the model's Fuller correlation at 298.15 K and 101,325 Pa. These are estimates, not necessarily measured binary diffusion coefficients.

<details>
<summary><strong>Show all 28 thermophysical and diffusivity records</strong></summary>


| Chemical | Formula | CAS No. | MW / (g/mol) | Cp / (J/mol/K) | Fuller V / (cm³/mol) | Dᵢ-air ×10⁵ / (m²/s) | Diffusion-volume provenance |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Hydrogen (H2) | H2 | 1333-74-0 | 2.0159 | 28.84 | 6.12 | 7.678 | tabulated Fuller volume for H2 |
| Nitrogen (N2) | N2 | 7727-37-9 | 28.0134 | 29.12 | 18.5 | 2.006 | tabulated Fuller volume for N2 |
| Oxygen (O2) | O2 | 7782-44-7 | 31.9988 | 29.38 | 16.3 | 2.023 | tabulated Fuller volume for O2 |
| Argon (Ar) | Ar | 7440-37-1 | 39.9480 | 20.79 | 16.2 | 1.929 | tabulated Fuller volume for Ar |
| Carbon monoxide (CO) | CO | 630-08-0 | 28.0101 | 29.14 | 18 | 2.024 | tabulated Fuller volume for CO |
| Carbon dioxide (CO2) | CO2 | 124-38-9 | 44.0095 | 37.13 | 26.9 | 1.595 | tabulated Fuller volume for CO2 |
| Methane (CH4) | CH4 | 74-82-8 | 16.0425 | 35.69 | 25.14 | 2.124 | estimated Fuller group contribution (C + H) |
| Ethane (C2H6) | C2H6 | 74-84-0 | 30.0690 | 52.51 | 45.66 | 1.430 | estimated Fuller group contribution (C + H) |
| Ethylene (C2H4) | C2H4 | 74-85-1 | 28.0532 | 42.85 | 41.04 | 1.515 | estimated Fuller group contribution (C + H) |
| Propane (C3H8) | C3H8 | 74-98-6 | 44.0956 | 73.60 | 66.18 | 1.137 | estimated Fuller group contribution (C + H) |
| n-Butane (n-C4H10) | C4H10 | 106-97-8 | 58.1222 | 98.49 | 86.7 | 0.969 | estimated Fuller group contribution (C + H) |
| Dimethyl ether (DME) | C2H6O | 115-10-6 | 46.0684 | 65.57 | 51.77 | 1.242 | estimated Fuller group contribution (C + H + O) |
| Ammonia (NH3) | NH3 | 7664-41-7 | 17.0305 | 35.65 | 20.7 | 2.227 | tabulated Fuller volume for NH3 |
| Hydrogen sulfide (H2S) | H2S | 7783-06-4 | 34.0810 | 34.20 | 27.52 | 1.671 | estimated Fuller group contribution (H + S) |
| Sulfur dioxide (SO2) | SO2 | 7446-09-5 | 64.0660 | 39.87 | 41.8 | 1.272 | tabulated Fuller volume for SO2 |
| Water vapor (H2O) | H2O | 7732-18-5 | 18.0153 | 33.58 | 13.1 | 2.537 | tabulated Fuller volume for H2O |
| Chlorine (Cl2) | Cl2 | 7782-50-5 | 70.9060 | 33.95 | 38.4 | 1.293 | tabulated Fuller volume for Cl2 |
| Isobutane (i-C4H10) | C4H10 | 75-28-5 | 58.1222 | 96.65 | 85.8 | 0.973 | Fuller atomic-increment estimate; not a measured binary diffusion volume |
| Ethylene oxide (C2H4O) | C2H4O | 75-21-8 | 44.0526 | 47.85 | 46.4 | 1.306 | Fuller atomic-increment estimate; not a measured binary diffusion volume |
| Hydrogen cyanide (HCN) | HCN | 74-90-8 | 27.0253 | 35.85 | 24.17 | 1.850 | Fuller atomic-increment estimate; not a measured binary diffusion volume |
| Nitrogen dioxide (NO2) | NO2 | 10102-44-0 | 46.0055 | 36.97 | 16.65 | 1.858 | Fuller atomic-increment estimate; not a measured binary diffusion volume |
| Hydrogen chloride (HCl) | HCl | 7647-01-0 | 36.4609 | 29.14 | 21.48 | 1.792 | Fuller atomic-increment estimate; not a measured binary diffusion volume |
| Carbon disulfide (CS2) | CS2 | 75-15-0 | 76.1410 | 45.66 | 50.5 | 1.154 | Fuller atomic-increment estimate; not a measured binary diffusion volume |
| 1,3-Butadiene (C4H6) | C4H6 | 106-99-0 | 54.0904 | 79.81 | 77.88 | 1.025 | Fuller atomic-increment estimate (C/H); not a measured binary diffusion volume |
| Methyl mercaptan (CH3SH) | CH4S | 74-93-1 | 48.1070 | 50.26 | 41.42 | 1.340 | Fuller atomic-increment estimate (C/H/S); not a measured binary diffusion volume |
| Phosphine (PH3) | PH3 | 7803-51-2 | 33.9976 | 37.11 | — | — | Not entered: no validated Fuller volume identified in the curated source set |
| Phosgene (COCl2) | CCl2O | 75-44-5 | 98.9160 | 57.68 | 60.98 | 1.038 | Fuller atomic-increment estimate (C/Cl/O); not a measured binary diffusion volume |
| Nitric oxide (NO) | NO | 10102-43-9 | 30.0061 | 29.86 | 11.17 | 2.313 | Fuller atomic-increment estimate (N/O); not a measured binary diffusion volume |


</details>

### 8.2 Flammability and toxicological endpoints

LFL and UFL values are in vol% in air. Toxicological endpoints are in ppm and retain their respective exposure bases. IDLH values based on 10% LEL are identified and excluded from toxicological contour screening.

<details>
<summary><strong>Show all 28 flammability and toxicity records (including database notes)</strong></summary>


| Chemical | LFL / (vol%) | UFL / (vol%) | AEGL-2, 60 min / (ppm) | ERPG-2, 1 h / (ppm) | IDLH / (ppm) | IDLH basis / interpretation | Database notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Hydrogen (H2) | 4 | 75 | — | — | — | — | No AEGL/ERPG/IDLH value was listed in the cited CAMEO entry. Hydrogen is an asphyxiant and a highly flammable gas. |
| Nitrogen (N2) | — | — | — | — | — | — | Simple asphyxiant; no chemical-specific acute toxic endpoint entered. |
| Oxygen (O2) | — | — | — | — | — | — | Oxygen is an oxidizer and can intensify combustion; it is not treated as a fuel component in Le Chatelier calculations. |
| Argon (Ar) | — | — | — | — | — | — | Simple asphyxiant; no chemical-specific acute toxic endpoint entered. |
| Carbon monoxide (CO) | 12.5 | 74 | 83 | 350 | 1200 | acute inhalation toxicity | AEGL-2 is for 60 min. ERPG-2 is a 1-hour guideline. IDLH from NIOSH. |
| Carbon dioxide (CO2) | — | — | — | — | 40000 | acute inhalation toxicity / asphyxiation | IDLH = 40,000 ppm. CO2 can cause hypercapnia and asphyxiation; no AEGL-2/ERPG-2 value entered here. |
| Methane (CH4) | 5 | 15 | — | — | 5000 | 10% LEL criterion; not a toxicological endpoint | NIOSH-listed IDLH is 10% LEL. Do not interpret this IDLH as an acute toxicity threshold. Methane is a simple asphyxiant. |
| Ethane (C2H6) | 3 | 12.5 | — | — | 3000 | 10% LEL criterion; not a toxicological endpoint | NIOSH-listed IDLH is 10% LEL. Do not interpret this IDLH as an acute toxicity threshold. Ethane is a simple asphyxiant. |
| Ethylene (C2H4) | 2.75 | 28.6 | — | — | — | — | No AEGL/ERPG/IDLH value entered from the cited CAMEO entry. |
| Propane (C3H8) | 2.1 | 9.5 | 17000 | — | 2100 | 10% LEL criterion; not a toxicological endpoint | EPA flags AEGL-2 (17,000 ppm) above 50% of its stated LEL; explosion hazard must be considered. NIOSH IDLH is based on 10% LEL, not toxicity. |
| n-Butane (n-C4H10) | 1.6 | 8.4 | 17000 | — | 1600 | 10% LEL criterion; not a toxicological endpoint | EPA flags AEGL-2 (17,000 ppm for 30 min to 8 h) above 50% of its stated LEL; explosion hazard must be considered. NIOSH IDLH is based on 10% LEL, not toxicity. |
| Dimethyl ether (DME) | 2 | 50 | — | — | — | — | CAMEO/USCG source reports LEL 2 vol% and UEL 50 vol%. Other compilations (for example NFPA-derived HSDB data surfaced through PubChem) report 3.4/27 vol%; this source discrepancy materially affects the contour. No AEGL/ERPG/IDLH value was available in the cited CAMEO entry. Select and document the source appropriate to the study. |
| Ammonia (NH3) | 15 | 28 | 160 | 150 | 300 | acute inhalation toxicity | AEGL-2 is for 60 min. ERPG-2 is a 1-hour guideline. NIOSH notes NH3 should be treated as flammable despite its DOT transport classification. |
| Hydrogen sulfide (H2S) | 4 | 44 | 27 | 30 | 100 | acute inhalation toxicity | AEGL-2 is for 60 min. ERPG-2 is a 1-hour guideline. NIOSH cautions that odor is not a reliable warning because olfactory fatigue occurs. |
| Sulfur dioxide (SO2) | — | — | 0.75 | 3 | 100 | acute inhalation toxicity | AEGL-2 is for 60 min. ERPG-2 is a 1-hour guideline. |
| Water vapor (H2O) | — | — | — | — | — | — | Water vapour is not treated as a fuel or as a toxic endpoint in this starter database. |
| Chlorine (Cl2) | — | — | 2 | 3 | 10 | acute inhalation toxicity | AEGL-2 is for 60 min. ERPG-2 is a 1-hour guideline. Chlorine is a strong oxidizer but is not itself a fuel gas. |
| Isobutane (i-C4H10) | 1.8 | 8.4 | — | — | — | — | CAMEO reports a flammable gas; no chemical-specific AEGL-2, ERPG-2 or IDLH acute-toxicity endpoint entered in this database. |
| Ethylene oxide (C2H4O) | 3 | 100 | 45 | 50 | 800 | NIOSH IDLH; potential occupational carcinogen noted by NIOSH | High consequence flammable and toxic gas. Verify current endpoint versions and chemical-specific handling before use. |
| Hydrogen cyanide (HCN) | 5.6 | 40 | 7.1 | 10 | 50 | acute inhalation toxicity | Highly toxic by inhalation. AEGL-2 is 60 min; ERPG-2 is 1 h; IDLH is a 30-min escape-planning reference. |
| Nitrogen dioxide (NO2) | — | — | 12 | 15 | 13 | acute inhalation toxicity | Toxic oxidizing gas; not a fuel. AEGL-2 is 60 min; ERPG-2 is 1 h. |
| Hydrogen chloride (HCl) | — | — | 22 | 20 | 50 | acute inhalation toxicity | Anhydrous hydrogen chloride gas. Do not confuse gas data with hydrochloric-acid solution data. |
| Carbon disulfide (CS2) | 1 | 50 | 160 | 50 | 500 | acute inhalation toxicity | Highly flammable and toxic vapour. The volatile-liquid source term is not represented by the current gas-only release module. |
| 1,3-Butadiene (C4H6) | 2 | 11.5 | 5300 | 500 | 2000 | 10% LEL criterion; not a toxicological endpoint | Interim AEGL-2 (60 min) is flagged by NOAA CAMEO as 10–49% of LEL; explosion hazard must also be considered. IDLH 2,000 ppm is based on 10% LEL and is excluded from toxic endpoint screening. Shipped as liquefied gas; polymerization and liquid release are not modeled here. |
| Methyl mercaptan (CH3SH) | 3.9 | 21.8 | 23 | 30 | 150 | acute inhalation toxicity | NOAA CAMEO identifies it as poisonous and flammable; odor is not a reliable warning. Low-boiling liquid/vapour behaviour and liquid spill evaporation are not modeled by the gas-only source term. |
| Phosphine (PH3) | 1.79 | — | 2 | 0.5 | 50 | acute inhalation toxicity | Highly toxic and pyrophoric. CAMEO reports LEL 1.79% but does not provide a numeric UEL; mixed-fuel Le Chatelier limits are therefore disabled when this species is present. No Fuller volume is supplied, so D in air is left missing rather than guessed. |
| Phosgene (COCl2) | — | — | 0.3 | 0.5 | 2 | acute inhalation toxicity | Extremely toxic choking agent. Gas or low-boiling volatile liquid; this database entry is for the gas-phase component and the existing source term does not model liquid phosgene evaporation or aerosol formation. |
| Nitric oxide (NO) | — | — | — | — | 100 | NIOSH IDLH; acute criteria and AEGL/ERPG not established in the cited sources | NIOSH lists IDLH 100 ppm and reports that NO is nonflammable but accelerates combustion. The cited NIOSH page does not provide LFL/UFL. Nitric oxide is rapidly converted in air to NO2; this simple single-species dispersion model does not model that chemistry. |


</details>

### 8.3 Database interpretation notes

- **DME:** the embedded LFL/UFL values are 2.0/50 vol% from the cited NOAA CAMEO/USCG source. Other compilations may report 3.4/27 vol%. This discrepancy materially affects flammability contours; select and document the source appropriate to the study.
- **Phosphine:** a numeric LFL is entered, but no numeric UFL is available in the cited record. The mixed-fuel Le Chatelier calculation is disabled if phosphine is included.
- **Nitric oxide:** treated as nonflammable but combustion-accelerating in the NIOSH source. It is not a fuel in Le Chatelier calculations. Its rapid conversion to NO₂ in air is not represented.
- **IDLH based on 10% LEL:** methane, ethane, propane, n-butane and 1,3-butadiene have IDLH values flagged as 10% LEL criteria. These are not toxicological endpoints and are excluded from toxicity contours.
- **Propane and n-butane:** the embedded EPA AEGL-2 values are flagged in the source notes as above 50% of their stated LEL. A toxic endpoint may therefore overlap a flammability concern; interpret both hazards together.
- **Diffusion volumes:** several organic-species volumes are Fuller group/atomic-increment estimates rather than measured values. PH₃ diffusivity is not entered because a defensible Fuller volume was not established in the curated source set.
- **Static properties:** MW, Cp and hazard data are static embedded records. Cp and thermophysical values are near ambient reference conditions. Verify current data, units, endpoint exposure duration, and local regulatory criteria before scenario use.

---

## 9. Assumptions, limitations, and validation

### 9.1 Assumptions and excluded physics

| Area | Current assumption / excluded effect |
|---|---|
| Gaussian dispersion | Passive gas; constant mean wind direction; horizontally homogeneous, time-invariant weather; Gaussian parameterisation; reflected ground; no terrain/building steering. |
| Meteorology | One Pasquill stability class and one roughness description throughout the domain; no evolving weather or wind direction. Very low wind/calm conditions are outside recommended use; the notebook warns when reference wind is below 1 m/s. |
| Source term | Single gas phase and ideal-gas equation of state; no pipe friction, vessel-wall heat transfer, jet trajectory, flashing, two-phase critical flow, droplet rainout or pool evaporation. |
| Mixture | Ideal-gas mixture; fixed composition; constant Cp near 298.15 K; no non-ideal interactions, chemical reaction, selective condensation or differential species loss. |
| Hazard screening | Le Chatelier is an approximate fuel-blend rule; no oxygen-enrichment/inerting correction; no additive toxic-mixture model; endpoint availability varies by component. |
| Dense gas | DEGADISx ground-area source only. Flat, level terrain and low source momentum. 3D display is profile post-processing, not native 3D solver output. |
| Numerics | Nonuniform grid and direct receptor comparison are diagnostics, not proof of grid convergence. Isosurface/contour levels are display settings unless explicitly matched to a documented endpoint. |

#### Specific source-term limitation for liquefied DME/LPG

The vessel blowdown model is for a rigid, adiabatic, well-mixed **single-phase ideal-gas** inventory. It does not include liquid phase, vapour-liquid equilibrium, real-gas behaviour, flashing, two-phase critical flow, aerosols, entrained droplets, composition fractionation, pipe friction, or pool evaporation. It should not be applied to a liquefied DME/LPG cylinder release unless the scenario is demonstrably a single-phase vapour event and the assumptions are justified.

#### Specific dense-gas limitation

The integrated DEGADISx pathway does not automatically calculate a high-pressure jet, pressurised-jet trajectory, flashing, droplets, aerosol rainout, pool spread, or evaporation before forming a ground-level dense-gas source. Do not feed a raw vessel/orifice discharge directly into the ground-area source model without a justified upstream source/jet/rainout model.

### 9.2 Recommended validation sequence

1. **Unit and implementation checks:** verify mixture-fraction conversion, wind-profile recovery at the reference height, release flow regime, monotonic blowdown pressure/temperature trends, finite/nonnegative concentration arrays, and source-profile mass balance.
2. **Grid convergence:** run coarse, medium, and fine meshes. Compare peak concentration, receptor histories, selected contour distances, and plume widths. Do not rely only on visual smoothness.
3. **Passive Gaussian validation:** compare with suitable tracer field data such as Prairie Grass using observed meteorology and stability classification. Compare centreline and crosswind profiles, downwind concentration, arrival time, and receptor concentration.
4. **Dense-gas validation:** run DEGADISx test cases or appropriate dense-gas field experiments. Validate the external DEGADISx backend separately from the notebook's 3D post-processor.
5. **Source-term validation:** compare ideal-gas orifice and blowdown outputs against independent compressible-flow calculations or measured single-phase vessel-discharge data. Do not validate a single-phase gas model against flashing/liquefied-gas tests without a two-phase source model.
6. **Hazard-data audit:** check each chemical property against current primary data. Record source edition/date, units, exposure duration, interim/final status, and IDLH basis.

Before applying the notebook to a scenario, document the selected model route, source term, weather basis, data provenance, grid-convergence outcome, benchmark comparison, and unresolved limitations. If a required physical process is not represented, add a validated sub-model or state that the scenario is outside the model's scope.

---

## 10. References

The following references underpin the equations, software modules, and data sources used in the model. Chemical-specific URLs are also embedded in the chemical database dictionaries in the notebook. “n.d.” indicates that the cited repository or database page does not provide a stable publication year.


1. Beals, G. A. (1971). *Guide to local diffusion of air pollutants*. Air Weather Service (MAC), U.S. Air Force.
2. Briggs, G. A. (1973). *Diffusion estimation for small emissions: Preliminary draft*. Atmospheric Turbulence and Diffusion Laboratory, Oak Ridge, Tennessee.
3. Britter, R. E. (1989). Atmospheric dispersion of dense gases. *Annual Review of Fluid Mechanics, 21*, 317–344.
4. Brutsaert, W. (1982). *Evaporation into the atmosphere: Theory, history, and applications*. D. Reidel.
5. Businger, J. A. (1973). Turbulence transfer in the atmospheric surface layer. In D. A. Haugen (Ed.), *Workshop on Micrometeorology*. American Meteorological Society.
6. Colenbrander, G. W. (1980). A mathematical model for the transient behavior of dense vapor clouds. Proceedings of the 3rd International Symposium on Loss Prevention and Safety Promotion in the Process Industries, Basel, Switzerland.
7. Colenbrander, G. W., & Puttock, J. S. (1983). Dense gas dispersion behavior: Experimental observations and model developments. International Symposium on Loss Prevention and Safety Promotion in the Process Industries, Harrogate, England.
8. Fuller, E. N., Schettler, P. D., & Giddings, J. C. (1966). New method for prediction of binary gas-phase diffusion coefficients. *Industrial & Engineering Chemistry, 58*(5), 18–27. https://doi.org/10.1021/ie50677a007
9. Hanna, S. R., Briggs, G. A., & Hosker, R. P. (1982). *Handbook on atmospheric diffusion*. U.S. Department of Energy, Technical Information Center.
10. Jones, R., Lehr, W., Simecek-Beatty, D., & Reynolds, R. M. (2013). *ALOHA® (Areal Locations of Hazardous Atmospheres) 5.4.4: Technical Documentation*. NOAA Technical Memorandum NOS OR&R 43. U.S. Department of Commerce / NOAA. https://repository.library.noaa.gov/view/noaa/2669
11. Lee, U. (2026). *degadisx* (Version 0.1.4) [Python software]. PyPI. https://pypi.org/project/degadisx/0.1.4/
12. Lee, U. (2026). *degadisx* [Source code]. GitHub. https://github.com/lyullee/degadisx
13. National Institute of Standards and Technology. (2025). *NIST Chemistry WebBook*, SRD 69. https://doi.org/10.18434/T4D303
14. National Institute of Standards and Technology. (n.d.). *Computational Chemistry Comparison and Benchmark Database (CCCBDB)*. https://cccbdb.nist.gov/
15. NOAA Office of Response and Restoration. (n.d.). *CAMEO Chemicals*. https://cameochemicals.noaa.gov/
16. Palazzi, E., De Faveri, M., Fumarola, G., & Ferraiolo, G. (1982). Diffusion from a steady source of short duration. *Atmospheric Environment, 16*(12), 2785–2790.
17. Pasquill, F. (1961). The estimation of the dispersion of windborne material. *Meteorological Magazine, 90*, 33–49.
18. Shapiro, A. H. (1953). *The dynamics and thermodynamics of compressible fluid flow* (Vol. I). Wiley.
19. Spicer, T. O., & Havens, J. A. (1989). *User’s guide for the DEGADIS 2.1 dense gas dispersion model* (EPA-450/4-89-019). U.S. Environmental Protection Agency. https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=2000J5GU.TXT
20. U.S. Environmental Protection Agency. (n.d.). *Air quality dispersion modeling: Alternative models (DEGADIS and SLAB)*. https://www.epa.gov/scram/air-quality-dispersion-modeling-alternative-models
21. U.S. Environmental Protection Agency. (n.d.). *About Acute Exposure Guideline Levels (AEGLs)*. https://www.epa.gov/aegl/about-acute-exposure-guideline-levels-aegls
22. U.S. Environmental Protection Agency. (n.d.). *Access Acute Exposure Guideline Levels (AEGLs) values*. https://www.epa.gov/aegl/access-acute-exposure-guideline-levels-aegls-values
23. VaibhavVasdev. (n.d.). *Gaussian-Plume_Model* [GitHub repository]. https://github.com/VaibhavVasdev/Gaussian-Plume_Model
24. National Institute for Occupational Safety and Health. (n.d.). *NIOSH Pocket Guide to Chemical Hazards*. Centers for Disease Control and Prevention. https://www.cdc.gov/niosh/npg/npg.html
25. COMSOL AB. (2024). *Transport properties: Gas diffusivity correlations* (COMSOL Multiphysics 6.4 documentation). https://doc.comsol.com/6.4/doc/com.comsol.help.lgp/lgp_ug_thermo.5.27.html


### Database source register

- **Thermophysical properties:** [NIST Chemistry WebBook (SRD 69)](https://doi.org/10.18434/T4D303) and [NIST CCCBDB](https://cccbdb.nist.gov/).
- **Flammability and chemical emergency data:** [NOAA CAMEO Chemicals](https://cameochemicals.noaa.gov/) and [NIOSH Pocket Guide to Chemical Hazards](https://www.cdc.gov/niosh/npg/npg.html).
- **AEGL endpoints:** [U.S. EPA AEGL programme](https://www.epa.gov/aegl/access-acute-exposure-guideline-levels-aegls-values).
- **ERPG values:** NOAA CAMEO records and referenced emergency-response guidance, as embedded in the database.
- **Dense-gas model:** [DEGADISx on PyPI](https://pypi.org/project/degadisx/0.1.4/), [DEGADISx source repository](https://github.com/lyullee/degadisx), and the [EPA alternative dispersion models page](https://www.epa.gov/scram/air-quality-dispersion-modeling-alternative-models).

## Repository files

A suggested repository layout is:

```text
.
├── README.md
└── 3D_Gas_Dispersion_Model.ipynb
```

## Disclaimer and contribution note

This README documents the current notebook implementation and its stated assumptions. Model outputs are estimates and should be independently checked before engineering use. Chemical-property and hazard-reference values can be revised by their source organisations; check the current primary record before use.
