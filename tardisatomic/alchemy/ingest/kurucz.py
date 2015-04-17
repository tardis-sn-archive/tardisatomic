import pandas as pd

from tardisatomic.alchemy import Atom, Ion
from tardisatomic.alchemy.ingest import BaseIngest

from tardisatomic.kurucz.io import (read_gfall_raw, parse_gfall,
                                    extract_levels, extract_lines)

class IngestGFAll(BaseIngest):


    def requirements_satisfied(self):
        if self.atomic_db.session.query(Ion).count() > 0:
            return True

    def ingest(self, fname):
        gfall = read_gfall_raw(fname)
        gfall_parsed = parse_gfall(gfall)
        levels = extract_levels(gfall_parsed)
        lines = extract_lines(gfall_parsed, levels)
