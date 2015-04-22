from tardisatomic.alchemy import Ion, Level, Transition, TransitionType, TransitionValue, DataSource, \
    TransitionValueType, Unit, DType
from tardisatomic.alchemy.ingest import BaseIngest

from tardisatomic.kurucz.io import (read_gfall_raw, parse_gfall,
                                    extract_levels, extract_lines)
from tardisatomic.util import convert_air_to_vacuum
import pandas as pd
import numpy as np


class IngestGFAll(BaseIngest):
    def requirements_satisfied(self):
        if self.atomic_db.session.query(Ion).count() > 0:
            return True

    def _convert_air_to_vacuum(self, wavelength):
        if wavelength > 2000:
            return convert_air_to_vacuum(self)
        else:
            return wavelength

    def ingest(self, fname):
        gfall = read_gfall_raw(fname)
        gfall_parsed = parse_gfall(gfall)
        levels = extract_levels(gfall_parsed)

        lines = extract_lines(gfall_parsed, levels)

        units_angstrom = Unit(unit='angstrom')
        units_none = Unit(unit='None')
        dtype_float = DType(short_name='float', name='float')
        data_source_ku = DataSource(short_name='KU', name='Kurucz_gfall')
        transition_type_line = TransitionType(name='Line', comment='')
        value_type_wl = TransitionValueType(name='wavelength', unit=units_angstrom, dtype_id=dtype_float)
        value_type_loggf = TransitionValueType(name='loggf', unit=units_none, dtype_id=dtype_float)
        self.atomic_db.session.add(units_angstrom)
        self.atomic_db.session.add(units_none)
        self.atomic_db.session.add(dtype_float)
        self.atomic_db.session.add(data_source_ku)
        self.atomic_db.session.add(transition_type_line)
        self.atomic_db.session.add(value_type_wl)
        self.atomic_db.session.add(value_type_loggf)

        levels['alchemy_level'] = 0
        _tmp = []
        for i, row in levels.iterrows():
            clevel = Level(ion_id=row['level_id'], level_number=row['level_number'], g=row['g'], energy=row['energy'],
                           label=row['label'], theoretical=row['theoretical'], data_source=data_source_ku)
            _tmp.append(clevel)
            row['alchemy_level'] = clevel

        levels['alchemy_leve'] = pd.Series(_tmp, index=levels.index)

        for i, row in lines.iterrows():
            lowerlevel = levels.ix[row['level_id_lower']]['alchemy_level']
            upperlevel = levels.ix[row['level_id_upper']]['alchemy_level']
            clineT = Transition(transition_type=transition_type_line, source_level=upperlevel,
                                target_level=lowerlevel, data_source=data_source_ku)
            clinegf = TransitionValue(transition=clineT, transition_value_type=value_type_loggf,
                                      value=row['loggf'], data_source=data_source_ku)
            vacum_wl = self._convert_air_to_vacuum(row['wavelength'])
            clinewl = TransitionValue(transition=clineT, transition_value_type=value_type_wl, value=vacum_wl,
                                      data_source=data_source_ku)
            self.atomic_db.session.add(lowerlevel)
            self.atomic_db.session.add(upperlevel)
            self.atomic_db.session.add(clineT)
            self.atomic_db.session.add(clinegf)
            self.atomic_db.session.add(clinewl)

        self.atomic_db.session.commit()


