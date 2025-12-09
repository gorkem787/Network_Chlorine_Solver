
# Chlorine Decay Simulation in Drinking Water Pipes

This project models the transport, dispersion and decay of chlorine within drinking water distribution systems using hydraulic data from an EPANET network. It numerically solves the advection-dispersion-reaction equations to simulate chlorine concentration over time and space.

## 🎯 Purpose

- Simulate the transport, dispersion and decay of chlorine and florine in pipes using EPANET hydraulic data.
- Analyze the effect of the decay coefficient (K) on chlorine concentration.
- Compare numerical results with EPANET's built-in chlorine simulation.

## 🧰 Requirements

- Python 3.7+
- Libraries:
  - `numpy`
  - `matplotlib`
  - `pandas`
  - `numba`
- EPANET Python Wrapper (`epyt`)
- `ADE_solver.py` – contains the numerical solution functions

Install required libraries:

```bash
pip install numpy matplotlib pandas epyt numba
```


## 🚀 Usage

To run the simulation, use the following command:

```bash
python ADE.py
```

A typical usage inside `ADE.py`:
Define pump as negative demand in epanet
```python
a = system(max_time=55 * 3600, time_step=3600 / 4)
nodes, pipes = a.constructor("epanet_net2/net2.inp", K=1.5 / 86400)
a.system_calculate_time()
a.plot_with_epanet(save=False)
```

## 📊 Outputs

- Time series plots of chlorine concentration in selected pipes.
- Comparison with EPANET's chlorine simulation results.
- Analysis of decay trends and system behavior.

## ✅ Features

- Read EPANET hydraulic and quality data.
- Divide pipes into customizable segments for numerical accuracy.
- Use first-order decay reaction kinetics.
- Simulate node mixing and boundary inflow/outflow.
- Generate detailed plots for analysis and comparison.
- For now, you can use with 1 reservoir and 1 pump.



---

Feel free to open an issue or submit a pull request if you’d like to contribute or have suggestions!
