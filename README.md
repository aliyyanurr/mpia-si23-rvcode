# MPIA Summer Internship 2023 - RV Code
Repository for RV code I created during my internship at MPIA, 2023. This code can be used for simple exoplanet and binary stars RV modeling.

### Libraries:
Before using the code, please install `kepler.py`:

```
pip install kepler.py
```

### Usage:
This code is better run in an interactive notebook `*.ipynb` file instead of running it in `*.py`.

Steps to use the code:
* Copy and run the code in `RV.py` to the first block of the Python notebook
* Define the input parameters: primary mass $M_1$ (kg), secondary mass $M_2$ (kg), semi-major axis $a$ (m), eccentricity $e$, inclination $i$ (rad), argument of periastron $\omega$ (rad), period $P$ (s), and the time of periastron passage $T_0$ (s). <br>
   You can convert $a$ to $P$ or vice versa using `ptoa` or `atop` function from `RV.py`.
* Create an array of time $t$ you want to simulate, with $T_0$ as the first element of the array. Ensure that all the time unit in the array is in seconds (s).
* Steps on simulation:
1. Calculate the _**average angular speed**_ of the (secondary) object in the orbit using `avg_ang_speed` function, using $P$ as the function's input.
2. Calculate the _**mean anomaly**_ using `mean_anomaly` function, using $t$ and _average angular speed_ as the input.
3. Calculate the _**eccentric anomaly**_ and _**true anomaly**_ using `kepler_solve` function, with _mean anomaly_ and $e$ as the input.
4. Calculate the polar coordinate position of the object in orbit, $r$ and $\theta$ with `pos` function, using $a$, $e$, and _true anomaly_ as the input.
5. Finally, calculate the radial velocity with `radial_velocity` function, using _true anomaly_ $\omega$, $e$, _systemic velocity_, $a$, $i$, and $M_1$, $M_2$.

Finally, You can plot RV vs $t$.

An example `ipynb` file for any RV simulation (2-body and multiple bodies) is provided in this repository.

Note: beware of the notation of the system when doing the simulation.

This code is a very simple form of simulating the RV plot and needs a lot of improvisation. Feel free to contribute to this repository!
