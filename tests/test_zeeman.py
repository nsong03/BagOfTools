import unittest
import numpy as np
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
from atomic_tools import zeeman


class TestZeeman(unittest.TestCase):
    def test_rb87_ground_gf(self):
        g_f_F1 = zeeman.rb87_ground_gf(1)
        g_f_F2 = zeeman.rb87_ground_gf(2)
        self.assertAlmostEqual(g_f_F1, -0.5018, places=4)
        self.assertAlmostEqual(g_f_F2, 0.4998, places=4)

    def test_energy_frequency_consistency(self):
        g_f = zeeman.rb87_ground_gf(2)
        B = np.linspace(0, 1e-3, 5)
        m_f = 1
        energy = zeeman.zeeman_shift(B, m_f, g_f)
        freq = zeeman.zeeman_frequency(B, m_f, g_f)
        np.testing.assert_allclose(energy / zeeman.HPLANCK, freq)


if __name__ == '__main__':
    unittest.main()
