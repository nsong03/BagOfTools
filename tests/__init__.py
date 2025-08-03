"""Test helpers for BagOfTools.

Ensures the repository root is on ``sys.path`` so local imports
like ``atomic_tools`` work regardless of the working directory.
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
