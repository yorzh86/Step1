<<<<<<< HEAD
Mark1
=====

The Mark 1 MD code is a developmental MD code with the intent of first duplicating the basic operation of Lammps.
After comparison with Lammps suggests that all basic operations are implemented correctly the real numbers will be replaced with dual numbers for automatic differentiation of the code's results.

Units
-----

The units commonly used for MD simulations vary, but are frequently the so-called metal units, as shown in the table below.

| Quantity          | Unit                        |
|-------------------|-----------------------------|
| mass              | grams/mole                  |
| distance          | Angstroms                   |
| time              | picoseconds                 |
| energy            | eV                          |
| velocity          | Angstroms/picosecond        |
| force             | eV/Angstrom                 |
| torque            | eV                          |
| temperature       | Kelvin                      |
| pressure          | bars                        |
| charge            | multiple of electron charge |
| dipole            | charge*Angstroms            |
| density           | gram/cm^dim                 |
| dynamic viscosity | Poise                       |
| electric field    | volts/Angstrom              |

The Mark 1 code currently does all calculations in standard SI units and relies on proper floating point arithmetic to handle the scale differences.
=======
Mark1
=====

The Mark 1 MD code is a developmental MD code with the intent of first duplicating the basic operation of Lammps.
After comparison with Lammps suggests that all basic operations are implemented correctly the real numbers will be replaced with dual numbers for automatic differentiation of the code's results.

Units
-----

The units commonly used for MD simulations vary, but are frequently the so-called metal units, as shown in the table below.

| Quantity          | Unit                        |
|-------------------|-----------------------------|
| mass              | grams/mole                  |
| distance          | Angstroms                   |
| time              | picoseconds                 |
| energy            | eV                          |
| velocity          | Angstroms/picosecond        |
| force             | eV/Angstrom                 |
| torque            | eV                          |
| temperature       | Kelvin                      |
| pressure          | bars                        |
| charge            | multiple of electron charge |
| dipole            | charge*Angstroms            |
| density           | gram/cm^dim                 |
| dynamic viscosity | Poise                       |
| electric field    | volts/Angstrom              |

The Mark 1 code currently does all calculations in standard SI units and relies on proper floating point arithmetic to handle the scale differences.
>>>>>>> 2abc27c8e73f839682eb484cedd80b29d229f06f
