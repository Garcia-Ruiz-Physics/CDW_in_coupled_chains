# CDW in Coupled Chains

MATLAB codes for simulating the total energy and phase diagrams of 1D coupled atomic chains, as presented in:

https://arxiv.org/abs/2602.23208  
*A. Garcia-Ruiz et al.*

---

## Overview

This repository contains MATLAB scripts used to generate the figures in the manuscript. The code computes total energies and phase diagrams of charge-density-wave (CDW) states in coupled 1D atomic chains as a function of:

- Rigidity
- Doping level
- Inter-chain coupling strength
- Number of chains
- Relative chain geometry (parallel / skewed)

---

## Contents

<details>
<summary><b>Figure 1 — Single chain</b></summary>

- Computes total energy of a single 1D atomic chain
- Produces CDW vs normal-state phase diagram
- Parameters: rigidity, doping

File:
- `Figure_1.m`

</details>

---

<details>
<summary><b>Figure 2 — Two parallel chains</b></summary>

- Total energy vs rigidity
- Phase diagrams for coupled chains:
  - rigidity vs inter-chain coupling
  - rigidity vs doping

File:
- `Figure_2.m`

</details>

---

<details>
<summary><b>Figure 3 — N parallel chains</b></summary>

- Phase diagrams for multiple coupled chains
- Includes:
  - rigidity vs inter-chain coupling
  - rigidity vs doping

File:
- `Figure_3.m`

</details>

---

<details>
<summary><b>Figure 4 — Two skewed chains</b></summary>

- Total energy for fixed rigidity (κ = 0.7)
- Phase diagrams:
  - rigidity vs inter-chain coupling
  - rigidity vs doping

File:
- `Figure_4.m`

</details>

---

<details>
<summary><b>Figure 5 — N skewed chains</b></summary>

- Phase diagram of skewed-coupled chains
- Parameters:
  - rigidity
  - doping

File:
- `Figure_5.m`

</details>

---

## Usage

Run each script in MATLAB:

```matlab
Figure_1
Figure_2
Figure_3
Figure_4
Figure_5
