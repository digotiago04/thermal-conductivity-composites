# 🔥 Effective Thermal Conductivity of Composites

This repository features MATLAB codes developed to compute the effective thermal conductivity of periodic composite materials. The implementations encompass the analytical formulation of the **Locally-Exact Homogenization Theory (LEHT)** and the numerical **Finite-Volume Theory (FVT)**, utilizing both the Mean-Field Theory approach and the energy-based approach.

The models are designed for composites consisting of an isotropic matrix and circular isotropic inclusions, represented through a periodic unit cell. In addition to determining the effective thermal conductivity matrix, the repository provides the visualization of the temperature field and the microscopic temperature profiles along the coordinate directions. Thus, this computational framework constitutes a robust tool for the comparative analysis between analytical and numerical homogenization methods. 

---

## 📐 Analytical solution

### Locally-Exact Homogenization Theory (LEHT)
The Locally-Exact Homogenization Theory (LEHT) is an analytical approach based on the Trefftz concept, in which the local fields are represented by series expansions that satisfy the governing differential equations. The solution is obtained by imposing continuity conditions at the fiber–matrix interface and periodicity conditions on the unit cell. This methodology allows the effective thermal conductivity of materials with inclusions to be determined. The LEHT formulation implemented in this repository is based on the concepts presented in **DOI:** [https://doi.org/10.1016/j.ijheatmasstransfer.2020.119477].

## NUMERICAL EXAMPLE
The example below shows how the LEHT.m script can be used. The input parameters are shown in Table 1.
* LEHT(k_m, k_i, frac, field, x_cut, y_cut)
**Table 1:** Inputs parameters' declaration
---
| Parameter | Type | Description |
| :---: | :---: | :--- |
| **`k_m`** | Float | Thermal conductivity of the matrix phase. |
| **`k_i`** | Float | Thermal conductivity of the inclusion phase. |
| **`frac`** | Float | Volume fraction of the inclusion in the Representative Unit Cell (RUC). |
| **`field`** | Integer (Flag) | Enables (`1`) or disables (`0`) the plotting of the total 2D temperature field. |
| **`x_cut`** | Float | Coordinate $x_1$ to extract the vertical temperature profile (cut). |
| **`y_cut`** | Float | Coordinate $x_2$ to extract the horizontal temperature profile (cut). |



---
## 🔲 Finite-Volume Theory (FVT)
FVT is a numerical approach based on the spatial discretization of the RUC into subvolumes (finite volumes). To calculate the effective thermal conductivity from the obtained local fields, this repository offers **two distinct mathematical formulations**:

* **Based on Mean-Field Theory:** In this formulation, homogenization is performed through the direct application of the volume averaging theorem. The effective conductivity is calculated by relating the volume average of the local heat fluxes to the average of the thermal gradients applied to the RUC.

* **Based on Energy Theory:** This formulation focuses on energy conservation. Homogenization is performed by ensuring the equivalence of the thermal energy dissipation rate (or entropy production rate) between the original heterogeneous RUC and the equivalent homogeneous material (effective medium).

---

##  💻 Requirements
The implementation of this tool was entirely developed in the MATLAB environment (version R2022b). Its development did not require the use of additional tools or packages, so the code can be executed in a standard MATLAB installation.

---

## 🚀 How to Use

All scripts were developed in a modular way in MATLAB, requiring no additional toolboxes to run the direct analyses.

### Input Parameters
At the beginning of each main script, the user can define the following geometric and material parameters:
* `L`, `H`: Dimensions of the Representative Unit Cell (RUC).
* `nx`, `ny`: Mesh discretization (for FVT) or number of terms in the expansions (for LEHT).
* `k_m`: Thermal conductivity of the matrix material.
* `k_i`: Thermal conductivity of the inclusion material.
* `frac`: Volume fraction of the inclusion in the composite.

### Execution
Simply open the desired script (e.g., `main_LEHT_isotropic.m` or the corresponding FVT files) in the MATLAB environment and press *Run*. 

### Results (Outputs)
The programs calculate and print directly to the Command Window the **Effective Thermal Conductivity Matrix ($K^*$)** of dimension $2 \times 2$, reflecting the macroscopic properties of the material in the $X$ and $Y$ directions.

---

## 📚 References and Theoretical Basis
The implemented mathematical models are based on advanced literature regarding computational micromechanics and periodic thermal homogenization. If you use these codes in your academic or research work, please consider citing the relevant bibliographic references associated with the development of these methods.
