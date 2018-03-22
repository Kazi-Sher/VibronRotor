# VibronRotor: Finite-Element Rotordynamic Code
Authors: Kazi Sher Ahmed (kazisherahmed@gmail.com) and [Prof. S. M. Ahmad]( https://www.giki.edu.pk/Faculty/Professor-S-M-Ahmad) (smahmad@giki.edu.pk)

License: VibronRotor is licensed under the terms of [GNU General Public License 3.0.]( https://www.gnu.org/licenses/gpl-3.0.txt)

## DESCRIPTION
VibronRotor is a finite-element code for prediction of the lateral rotordynamic response of flexible rotors. Employed finite-element formulation is based on the work of Nelson and McVaugh [1]. This code provides analysis tools which enable an appropriate selection of rotor design parameters for stable operation. Mesh approach in the code limits the element length-to-diameter ratio within the user-provided value for modeling accuracy. Response prediction, on a basic level, relies on eigenanalysis and steady-state imbalance response analysis. 

## USAGE
A free scientific programming language [GNU Octave](https://www.gnu.org/software/octave/#install) is the preferred environment to run the code. However, due to bidirectional syntactic compatibility of GNU Octave with MATLAB, users can also execute the code on MATLAB.

*core.m* is the parameter input file to define rotor geometrical / mechanical properties, speed-dependent bearing coefficients, and functionalities controls. Selected analysis is executed once *core.m* is run. An interactive GUI is in development to replace the parameter input method.

## FUNCTIONALITIES WITH ROTOR DESIGN INSIGHTS

| Functionalities | Design Insights [2-4] |
| ------------------ |:----------:| 
| Mode shapes | <p align="left"> When plotted for a range of bearing stiffness values, mode shapes depict the rigidity or flexibility of rotor during critical speed excitation. Mode shapes also guide the placement of imbalance masses in response analysis in accordance with API and ISO standards. |
| Campbell diagram (CD) | <p align="left"> In CD, intersection of bifurcated damped natural frequencies with synchronous excitation line reveals the critical speeds of a rotor system. CD assessment, in conjunction with modification in bearing stiffness, spans, and rotor component mass properties, guides the placement of critical speeds relative to rotor operating speed. Additionally, comparison of CD with imbalance response plots identifies the sensitivity of modes to residual imbalance. |
| Critical speed map (CSM) | <p align="left"> CSM guides the change in bearing stiffness values to ensure a suitable safety margin between critical speeds and system operating speed. |
| Imbalance response amplitude and phase | <p align="left"> Severity of rotor critical speed vibrations and hence the safety of system, against a range of bearing damping values, is gauged by imbalance response plots. |
| Orbit plot | <p align="left"> Construction of whirl orbits at various rotor speeds and rotor axial locations reveal the forward or backward directions of rotor whirl. Further, major axis of orbit ellipses can be compared with available rotor-stator clearance. |
| Instability threshold map | <p align="left"> Instability analysis reveals the lowest rotor speed and power where rotor is susceptible to self-excited vibrations. Machine operation around these speeds is avoided. |
| Mesh plot | <p align="left"> Here,  element division is imposed on the rotor schematic. | 

## REFERENCES
[1] Nelson H and McVaugh J. The dynamics of rotor-bearing systems using finite elements. Journal of Engineering for Industry 1976; 98: 593-600.

[2] Adams ML. Rotating machinery vibration: from analysis to troubleshooting. CRC Press, 2009.

[3] Chen WJ and Gunter EJ. Dynamics of rotor-bearing systems. Trafford publishing, Canada, 2010.

[4] Vance JM, Zeidan FY and Murphy B. Machinery vibration and rotordynamics. John Wiley & Sons, 2010.
