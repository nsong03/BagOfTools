# BagOfTools

A collection of utilities for atomic physics simulations. The toolkit is
in an early stage and currently contains code for computing Zeeman energy
splittings in Rubidium.

## Usage

The main functionality lives in `atomic_tools.zeeman`. Example:

```python
from atomic_tools import zeeman

# g_F for $^{87}$Rb ground state F=2
gF = zeeman.rb87_ground_gf(2)

# Frequency shift for m_F=2 in a 1 gauss magnetic field
B = 1e-4  # Tesla
mf = 2
shift_hz = zeeman.zeeman_frequency(B, mf, gF)
print(shift_hz)
```

Run tests with:

```bash
python -m unittest
```
