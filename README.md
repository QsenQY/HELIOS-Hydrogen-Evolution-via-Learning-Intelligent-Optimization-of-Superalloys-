# HELIOS Hydrogen-Evolution-via-Learning-Intelligent-Optimization-of-Superalloys
Machine learning drives automated modeling of high-entropy alloys and prediction of hydrogen production performance

![Workflow Diagram](workflow.png)

## Abstract
High-entropy alloy materials demonstrate exceptional catalytic properties due to their distinctive multi-component attributes and electronic effects. Nonetheless, the extensive data landscape of high-entropy alloys presents substantial hurdles in identifying high-performance catalysts. Compounding this challenge is the anisotropic character of surface sites within each catalyst, rendering performance prediction exceedingly difficult with conventional computational techniques. Although contemporary machine learning models have achieved significant advances in predicting properties for specific principal element combinations of high-entropy alloys, they encounter difficulties when modeling random combinations of principal elements and their concentrations. In this context, we introduce a closed-loop high-throughput workflow that integrates high-throughput automated modeling and performance prediction of high-entropy alloys and adsorption structures, along with a high-throughput synthesis method leveraging microchannel technology. We apply a previously established thermodynamic model to forecast the stability of principal element combinations and harness the MatterGen high-throughput modeling approach to screen 133,233 thermodynamically stable high-entropy alloy materials. Utilizing the high-throughput automated workflow, we have generated over 1,160,000 hydrogen atom adsorption models and, for the first time, have successfully predicted hydrogen adsorption energies using the EquiformerV2 model, pinpointing five high-entropy alloy materials with superior catalytic performance. These materials have been successfully synthesized via microchannel technology, and electrochemical experiments have confirmed their outstanding catalytic properties (overpotential = 5.5–9 mV). Statistical analysis indicates that the performance of active sites on high-entropy alloy surfaces adheres to a bimodal distribution. The unique electronic and structural synergy in high-entropy alloys substantially diminishes the inherent properties of principal elements, leading to localized averaging effects. This research presents innovative concepts and methodologies for navigating the data space of high-entropy alloy materials and for high-throughput prediction of catalytic performance.


## Download Data


- **v1.1** Release:  
  - [Data.zip (all models and adsorption‐energy statistics)](https://github.com/your-username/your-repo/releases/download/v1.1/Data.zip)

## External Dependencies

This project uses Git submodules to integrate two external repositories. Before proceeding, clone and initialize all submodules:

```bash
git clone --recurse-submodules git@github.com:QsenQY/HELIOS-Hydrogen-Evolution-via-Learning-Intelligent-Optimization-of-Superalloys-.git
cd HELIOS-Hydrogen-Evolution-via-Learning-Intelligent-Optimization-of-Superalloys-
git submodule update --init --recursive
```

### fairchem

- **Documentation** (upstream):  
  https://github.com/facebookresearch/fairchem#readme

- **Source (pinned)**:  
  https://github.com/QsenQY/fairchem/tree/977a80328f2be44649b414a9907a1d6ef2f81e95


### mattergen

- **Documentation** (upstream):  
  https://github.com/microsoft/mattergen#readme

- **Source (pinned)**:  
  https://github.com/QsenQY/mattergen/tree/ec029d177c93709fa9a2ea4e48b872760d09c63b
