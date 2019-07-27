# NMR utils

This python module implements and integrates equation number 5 from: Allard, P., Helgstand, M. and Hard, T.,
Journal of Magnetic Resonance, 129:19-29 (1997) for simulation of the mixing part of ROESY experiments.
It utilizes the Boost 1.68.0 odeint library for the integration and of course the Python/numpy C-APIs.
Once built with the "python3 setup.py install" command it will form a resident module of the python interpreter.
So far it provides one function and one function only, the "allard97_roesy_mixer" this function describe the 
cross relaxation between two nuclear spins during ROESY NMR spin-locking experiments. The "allard97_roesy_mixer"
function takes no less than 16 arguments (all floats), like this:

allard97_roesy_mixer(mixt, offset_a, R1a, R2a, offset_b, R1b, R2b, sl_strength, sigma_l, simga_t, MAx0, MAy0,
MAz0, MBx0, MBy0, MBz0)

Where:
    mixt = mixing time/spin-lock duaration (sec)
    offset_a = offset spin A relative to spin-lock carrier (Hz)
    R1a = longitudinal relaxation rate constant spin A (Hz)
    R2a = transeverse relaxation rate constant spin A (Hz)
    offset_b = offset spin B relative to spin-lock carrier (Hz)
    R1b = longitudinal relaxation rate constant spin B (Hz)
    R2b = transeverse relaxation rate constant spin B (Hz)
    sl_strength = spin-lock field strenght in (Hz)
    sigma_l = longitudinal cross-relaxation constant between spin A and B (Hz)
    sigma_t = transverse cross relaxation constant spin A and B (Hz)
    MAx0 = x-mangetization of spin A at at start of spin-lock
    MAy0 = y-mangetization of spin A at at start of spin-lock
    MAz0 = z-mangetization of spin A at at start of spin-lock
    MBx0 = x-mangetization of spin A at at start of spin-lock
    MAy0 = y-mangetization of spin A at at start of spin-lock
    MAz0 = z-mangetization of spin A at at start of spin-lock
 
It returns a numpy 2D array with 7 columns (time,MAx,MAy,MAz,MBx,MBy,MBz) containing the evolution of
magnetization during the spin-lock.


WARNING AND INFORMATION: This is a hobby project written by a non skilled programmer just for fun
and for the joy of programming. The code does compile and executes with reasonable result and for the moment
I don't think there are any memory leaks.

## Getting Started

The NmrModule.cpp file gets compiled with setup.py file found in this directory.

### Prerequisites

The c++ part utilize the Boost odeint 1.68.0 library so that is a necessity.
