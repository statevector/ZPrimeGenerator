# Z' Event Generator

This program is a physics event generator.
It computes the e+ e- --> mu+ mu- scattering cross section at the Leading-Order (LO) using Monte Carlo integration.
Events are generated according to the specific MC sampling of the cross section, which defines the kinematics of the final-state muons. 
This allows a large sample of events to be generated, which can then be used for subsequent statistical analysis (e.g. histogramming, etc.).

It is well established that the scattering process above can proceed through the exchange of the photon or Z boson of the Standard Model. 
But are these the only options nature has chosen? Many theories Beyond the Standard Model (BSM) predict at least one additional heavy, neutral gauge boson, known as the Z' boson, in complete analogy with the Standard Model Z boson.
If such a particle exists, then the above scattering process can proceed through Z' exchange as well!

Indeed, this process was searched for at the LEP collider and will be searched for again at upgraded future electron-positron colliders, such as CLIC or the ILC if they are built.
But how can such a particle be detected? The presence of a new, heavy Z' boson can be infered from subtle deviations in the angular distributions of the final-state muons from their Standard Model expectation (i.e. if only the photon and Z boson participate in the interaction)
These deviations are quantified in an observable called the forward-backward asymmetry (A_{FB}), and are ultimately due to the parity violating effects of the weak interaction.

The event generate here can be used to study these effects in light of a new Z' boson, whose specific properties, such as mass, decay width, and fermionic couplings can be set by the user.
These properties characterize the A_{FB} distribution and can even be used to distinguish one Z' model from another.
Therefore, the smoking-gun signature of Z' production would appear as a deviation in the A_{FB} distribution, compared to what the Standard Model predicts.
Detailed SM/BSM comparisons can be made with this code as running with a Z' is not required.

## How to Compile

g++ -o main.exe main.cpp  mt19937-64.c

## How to Run

From the command line, 
```
main.exe [\sqrt{s}]
```
where,
- \sqrt{s} is the center-of-mass energy of the electron-positron collider

## Examples

To run a simulation with \sqrt{s} = 500 GeV
```
main.exe 500
```
