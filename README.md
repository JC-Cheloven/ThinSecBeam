# ThinSecBeam
Thin-walled Sections of Beams

Certainly the crown jewel of this repository. ThinSecBeam offers functionality which is difficult to find even in top-level proprietary software. It analyzes the cross-section of a thin-walled beam under linear elastic hypothesis using the common mathematical model for this kind of beam. Closed paths in the section are supported and are in fact irrelevant in the adopted implementation.

It calculates nearly everything conceivable for a thin-walled section in the linear elastic domain, including warping modulus, torsion stiffness, the core, shear centre and other properties of the section. Additionally, if any of the stress resultants (shear force, bending moment, axial force) are provided, some plots of stress fluxes and warping displacements will be drawn. Among them, 3D practicable plots of warping displacements are of high value in understanding the deformation of the beam.

A hardcoded example and some example files are included.

ThinSecBeam is programmed in Python 3 and uses NumPy & Matplotlib as principal dependencies. Currently (as for Dec 2021) the Spanish and English versions of
the interface are available but the english version of the manual isn't ready yet. Coming soon.

Enjoy !!
