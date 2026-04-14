# MMDM
A simple, precise magnetic force formula was proposed by modifying the Magnetic Dipole Method with two correction parameters α and β. It applies to cylindrical and cuboidal magnets in nonlinear piezoelectric energy harvesters, achieving &lt;1% relative error even at small separations while maintaining analytical simplicity for theoretical studies.
# Precise magnetic force model for nonlinear piezoelectric energy harvesting

This repository provides MATLAB code that implements the **highly accurate magnetic force model** proposed in:

> Yang, Y., Xiang, H. (2023). **A simple and precise formula for magnetic forces in nonlinear piezoelectric energy harvesting**. *Nonlinear Dynamics*, 111(7), 6085–6110.  
> https://doi.org/10.1007/s11071-022-08160-5

The code calculates the voltage frequency response of a monostable piezoelectric energy harvester with magnetic end coupling. It compares the traditional magnetic dipole model (MDM) with the **present model**, demonstrating that the proposed model yields much more accurate results (relative error < 1% in typical cases).

## Results

- **Magnetic force curves** – the corrected model matches theoretical/numerical integration results closely, while the MDM shows significant deviation.
- **Voltage frequency response** – the corrected model accurately predicts the peak voltage and the overall shape of the response, in good agreement with finite element simulations (e.g., ANSYS) and experimental data.

## Run the code

1. Open MATLAB (R2023b or later, with Symbolic Math Toolbox).
2. Run `main.mlx`.
3. The script will generate two key figures:
   - Magnetic force comparison (MDM vs present model vs reference data)
   ![Magnetic force comparison](Figures/Magnetic_force.png)
   - Voltage frequency response (stable/unstable branches, comparison between MDM and present model)
   ![Voltage frequency response](Figures/FRC.png)

Optional Excel files (`Magnetic force_Upadrashta.xlsx` and `FR_Upadrashta.xlsx`) can be placed in the same folder to plot reference data points.

## Magnetic force correction functions

This repository provides two standalone MATLAB functions to compute the correction parameters α and β for cylindrical and cuboidal magnets:

- `correct_cylindermagnet_coeff(MA, R, T, d)` – for cylindrical magnets  
- `correct_cuboidalmagnet_coeff(MA, W, H, L, d)` – for cuboidal magnets (square cross-section recommended)

**Inputs:**  
- `MA` : magnetization (A/m)  
- `R` / `W, H, L` : magnet half‑dimensions (m)  
- `T` : half‑thickness (m)  
- `d` : spacing between magnet centers (m) – can be a **vector**; the outputs will be arrays of the same length.

**Outputs:**  
- `alpha`, `beta` : correction coefficients (array if `d` is a vector)  
- `dmin` : minimum valid spacing (m)

**Performance:**  
The functions use **parallel computing (`parfor`)** to speed up calculations when `d` is a vector. Ensure the Parallel Computing Toolbox is installed.

These functions can be used independently in any MATLAB script:

```matlab
d_vec = 10e-3:1e-3:20e-3;
[alpha, beta, dmin] = correct_cylindermagnet_coeff(891e3, 5e-3, 2.5e-3, d_vec);
F_corrected = alpha .* beta.^2 .* F_MDM(beta .* w_tip);
```

See the main script `main.mlx` for a complete example of applying them to a nonlinear energy harvester.

## Citation

If you use this code, please cite the paper:

```bibtex
@article{yang2023simple,
  title={A simple and precise formula for magnetic forces in nonlinear piezoelectric energy harvesting},
  author={Yang, Yi and Xiang, Hongjun},
  journal={Nonlinear Dynamics},
  volume={111},
  number={7},
  pages={6085--6110},
  year={2023},
  publisher={Springer},
  doi={10.1007/s11071-022-08160-5}
}
