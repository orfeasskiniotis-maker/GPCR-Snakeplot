"""Parse GPCRdb generic residue numbers.

GPCRdb shows numbers like "3.50x50": the part before "x" is helix.Ballesteros
(class-specific), the part after "x" is the structure-based GPCRdb number we
use as the key. So "3.50x50" -> key "3x50", helix 3, position 50.

Bulges add a third digit: "2.55x551" means an extra residue between x55 and
x56, so its position is 55.1. Constrictions simply skip a number, so they
need no special handling: that slot stays empty.
"""
from typing import NamedTuple, Optional


class Generic(NamedTuple):
    key: str          # "3x50", the template key
    helix: int        # 3  (8 = helix 8; loops use two digits, e.g. 12 = ICL1)
    position: float   # 50.0, or 55.1 for a bulge


def parse(display: Optional[str]) -> Optional[Generic]:
    """'3.50x50' -> Generic('3x50', 3, 50.0); None for residues without a number."""
    if not display or "x" not in display:
        return None
    left, pos = display.split("x")
    helix = int(left.split(".")[0])
    # Two digits are the position; a third digit marks a bulge (551 -> 55.1).
    position = int(pos[:2]) + (int(pos[2:]) / 10 if len(pos) > 2 else 0.0)
    return Generic(f"{helix}x{pos}", helix, position)


if __name__ == "__main__":
    for s in ["3.50x50", "2.55x551", "45.50x50", None]:
        print(s, "->", parse(s))
