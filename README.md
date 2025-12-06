# 🍦 POpSiCLE

**PhOtoelectron SpeCtrum library for Laser-matter intEractions**

Welcome to **POpSiCLE**! 👋

POpSiCLE is a parallel library designed to facilitate the calculation of photoelectron spectra for various grid-based solutions of the time-dependent Schrödinger equation and the time-dependent Kohn-Sham equations of TDDFT.

## 🚀 Features & Methods

This library implements a variety of methods to calculate photoelectron spectra, giving you flexibility depending on your simulation needs. The methods implemented include:

*   **Time-Dependent Surface Flux Method (t-SURFF)**
*   **Sampling Point Method**
*   **Fourier Transform Methods** (transforming the spatial wavefunction into momentum space)
*   **Hankel Transforms**

## 🔌 Portability & Integration

POpSiCLE has been written to be as portable as possible! 🌍

It can be interfaced straightforwardly with other codes used for modelling laser-matter interactions, including:
*   Codes solving the Kohn-Sham equations of Time-Dependent Density Functional Theory (TDDFT).
*   Codes developed for specific research (such as thesis work).
*   Various other highly-parallelised codes.

## 📘 Theoretical Basis

### Calculation of Photoelectron Spectrum
Photoelectron spectroscopy is a powerful tool used to obtain valuable information about matter irradiated by intense laser pulses. By measuring the energy and angular distributions of the photoproducts emitted during interaction with these laser pulses, we can obtain valuable information about the structure and dynamical response of the material.

Calculating photoelectron spectra is computationally demanding. **POpSiCLE** implements three main approaches to solve this:

1.  **The Fourier Method**: The Fourier transform of the spatial scattering solution into momentum space.
2.  **The Sampling Point Method**: Collecting outgoing wavepackets at detector points.
3.  **The Time-Dependent Surface Flux Method (t-SURFF)**: Recording flux through a surface.

#### 1. Fourier Method
The simplest way to extract asymptotic information is by projecting the wavefunction onto field-free plane waves at the end of a calculation. This "brute force" method requires an interaction volume large enough to hold high-energy electrons. The desired spectra are obtained by performing a Fourier transform of the unbound wavefunction at the end of the pulse.

#### 2. The Sampling Point Method
To reduce computational demand, this method mimics experimental detection. We collect outgoing electron wavepackets as they pass certain "detector points" in coordinate space and energy-analyse this information. The calculation size is smaller since absorbing boundaries can be used to remove ionised wavepackets that have passed the detector points.

#### 3. The Surface Flux Method (t-SURFF)
The t-SURFF method records the flux of ionised electrons passing a surface of radius $r_b$. The spectra are obtained by calculating the projection of the outgoing electron wavepackets at this surface onto Volkov waves (solutions of a free electron in a laser field). This converts the volume integral (used in the Fourier method) into a time integral of a surface flux, dramatically reducing computational cost.

## ⚙️ Numerical Implementation

**POpSiCLE** is written in **FORTRAN 2003** and uses **MPI** for parallel communication and **HDF5** for efficient I/O. It is designed to interface with various coordinate systems (Cartesian, cylindrical, spherical).

### Fourier Transforms
The library provides efficient parallel routines for:
*   **Fourier Transform**: Using FFT (Fast Fourier Transform) for Cartesian coordinates.
*   **Hankel Transform**: For cylindrical coordinates (Fourier-Bessel transform).
It supports both **Intel MKL** and **FFTW** libraries.

### Interpolation
For surface-based methods (t-SURFF, Sampling Point), the wavefunction is interpolated onto spherical shells.
*   **Bi/Tri-cubic Interpolation**: The default method. Fast and accurate using piecewise polynomials.
*   **Shepard’s Method**: Available for scattered data, though significantly slower.

### Included Tools
*   `tsurff_calculator`: Calculates spectral amplitudes from surface files using the t-SURFF method.
*   `sampling_calculator`: Calculates PES using the sampling point method.
*   `popsicle_viewer`: A Python script for visualizing the results.

## 📜 Background & Credits

**POpSiCLE** was co-designed and developed to provide a robust solution for the physics community.

*   **Funding:** This work was performed as a job within the **eCSE program** funded by **ARCHER** (the UK National Supercomputing Service).
*   **Development:** Co-designed and developed as a library to allow the calculation of photoelectron spectra for various grid-based solutions.

---
*Happy coding!* ⚛️
