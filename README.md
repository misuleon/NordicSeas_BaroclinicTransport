# NordicSeas_BaroclinicTransport

This repository contains MATLAB scripts to compute the annual baroclinic transport into the Nordic Seas (0–500 m) using hydrographic profiles from the Rockall Trough and Norwegian Sea (1900–2023).

## Contents

- `compute_baroclinic_transport_NS_RT.m` – main script for transport time series
- `deseason.m` – function to remove seasonal cycle
- `sum_Z.m` – trapezoidal vertical integral
- `inside3.m` – polygon mask for spatial filtering
- `README.md` – this file

## How to Run

1. Load profiles (`dyn_RT`, `dyn_NS`) and depths (`dpth`)
2. Run `compute_baroclinic_transport_NS_RT.m`
3. Output: annual transport time series (`Trnsprt`) and uncertainties (`Trns_uncert`)

## Citation

If you use this code, please cite:

Chafik, L. (2025), *Baroclinic transport into the Nordic Seas from historical hydrography (1900–2022)*, GitHub v1.0.0, https://github.com/misuleon/NordicSeas_BaroclinicTransport

and cite the associated papers:

Chafik et al. (2025): **TBA**

Rossby et al. (2020):  https://doi.org/10.1029/2020GL087456

## License

MIT License
