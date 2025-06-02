# Nonplanar TMS Coil Placement SimNIBS

A numerical optimization toolkit for transcranial magnetic stimulation (TMS) coil placement on nonconvex target surfaces using double-cone coil geometries.

## 📄 Citation

This repository contains the implementation for the following research:

> Zhang, X., Hancock, R., Santaniello, S. (2025). Projection-Based Numerical Optimization of Transcranial Magnetic Coil Placement for Nonconvex Target Surfaces and Double-Cone Coil Geometries. *Journal of Neural Engineering*.

## 🔧 Dependencies

- **SimNIBS 4.1** - Simulation of non-invasive brain stimulation
- **MATLAB** - For running optimization scripts
- **SPM12** - Statistical parametric mapping

## 🚀 Getting Started

### Preparation

#### 1. Atlas Setup
Download the required atlases (T1 and T2 images) from the [COBRA Lab](https://www.cobralab.ca/atlases).

**Important:** Convert atlas files to NIFTI format (*.nii) using [MOMinc](https://github.com/SIMEXP/mominc) before proceeding.

Once converted, place the image files in the `raw/` directory.

#### 2. SimNIBS Configuration
Copy the following files to your SimNIBS installation:
- `opt_struct.py`
- `ADMlib.py`

**Destination:** `SIMNIBSDIR/optimization/`
- Where `SIMNIBSDIR = <SimNIBS installation path>/simnibs_env/Lib/site-packages/simnibs`

⚠️ **Backup the original files before replacing them.**

### Execution Workflow

#### Step 1: Prepare Atlas Data
```bash
# Run in MATLAB
add_neck_to_atlases.m
```

Then execute the CHARM segmentation pipeline in SimNIBS:
```bash
charm brain1_charm brain1_t1_extended.nii brain1_t2_extended.nii --registerT2
```

#### Step 2: Generate Coil Geometries
```matlab
% Run in MATLAB
get_coilpos_base.m
```
This generates resampled coil geometries for optimization.

#### Step 3: Generate Candidate Placements
```matlab
% Run in MATLAB
tms_optimization_presampling.m
```

**Configuration:** Modify parameters in the `%% Parameters` section to customize the optimization process.

#### Step 4: Optimize Coil Placement
```matlab
% Run in MATLAB
tms_optimization_run.m
```

**Options:** Choose between ADM (Alternating Direction Method) or direct optimization methods.

## 💡 Key Features

- **Nonplanar surface optimization** - Handles complex brain geometries
- **Double-cone coil support** - Specialized for advanced coil designs
- **Multiple optimization methods** - ADM and direct approaches
- **Cerebellum targeting** - Supports regions not covered by SimNIBS 4.5

## 📋 Notes

While SimNIBS 4.5 has implemented TMS coil placement optimization for nonplanar coils, it does not yet support certain brain regions such as the cerebellum (for lack of gifti surface data from CHARM segmentation). This toolkit provides a valuable alternative for such specialized applications.

## 🤝 Contributing

For questions, issues, or contributions, please refer to the original research paper or contact the authors.

## 📂 Repository Structure

```
├── raw/                    # Atlas image files (T1, T2)
├── opt_struct.py          # SimNIBS optimization structure
├── ADMlib.py             # ADM optimization library
├── add_neck_to_atlases.m # Atlas preparation script
├── get_coilpos_base.m    # Coil geometry generation
├── tms_optimization_presampling.m  # Candidate placement generation
├── tms_optimization_run.m          # Main optimization script
└── # Other MATLAB files
```