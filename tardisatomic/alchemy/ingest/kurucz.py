from tardisatomic.alchemy import Ion, Level, Transition, TransitionType, TransitionValue, DataSource, \
    TransitionValueType, Unit, DType, Ion, Atom
from tardisatomic.alchemy.ingest import BaseIngest

from tardisatomic.kurucz.io import (read_gfall_raw, parse_gfall,
                                    extract_levels, extract_lines)
from tardisatomic.util import convert_air_to_vacuum
import pandas as pd

import ipdb


class IngestGFAll(BaseIngest):
    def requirements_satisfied(self):
        if self.atomic_db.session.query(Ion).count() > 0:
            return True

    def _convert_air_to_vacuum(self, wavelength):
        if wavelength > 2000:
            return convert_air_to_vacuum(wavelength)
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
        value_type_wl = TransitionValueType(name='wavelength',
                                            unit=units_angstrom, dtype=dtype_float)
        value_type_loggf = TransitionValueType(name='loggf',
                                               unit=units_none, dtype=dtype_float)
        self.atomic_db.session.add(units_angstrom)
        self.atomic_db.session.add(units_none)
        self.atomic_db.session.add(dtype_float)
        self.atomic_db.session.add(data_source_ku)
        self.atomic_db.session.add(transition_type_line)
        self.atomic_db.session.add(value_type_wl)
        self.atomic_db.session.add(value_type_loggf)
        self.atomic_db.session.commit()

        levels['alchemy_level'] = 0
        _tmp = []
        for i, row in levels.iterrows():
            ion = self._find_ion(row['atomic_number'], row['ion_number'])
            clevel = Level(ion = ion,
                           level_number=row['level_number'], g=row['g'],
                           energy=row['energy'], label=row['label'],
                           theoretical=row['theoretical'],
                           data_source=data_source_ku)
            _tmp.append(clevel)
            row['alchemy_level'] = clevel

        levels['alchemy_level'] = pd.Series(_tmp, index=levels.index)
        levels['atomic_number'] = levels['atomic_number'].astype(int)
        levels['ion_number'] = levels['ion_number'].astype(int)
        levels['level_number'] = levels['level_number'].astype(int)
        indexed_levels = levels.set_index(['atomic_number', 'ion_number', 'level_number'])

        for i, row in lines.iterrows():
            catom_number = int(row['atomic_number'])
            cion_number = int(row['ion_number'])
            level_number_lower = int(row['level_number_lower'])
            level_number_upper = int(row['level_number_upper'])
            #print(
            #'atom: {0}, ion: {1}, level_l: {2}, level_u: {3})'.format(catom_number, cion_number, level_number_lower,
            #                                                          level_number_upper))
            lowerlevel = indexed_levels.loc[catom_number, cion_number, level_number_lower]['alchemy_level'].values[0]
            upperlevel = indexed_levels.loc[catom_number, cion_number, level_number_upper]['alchemy_level'].values[0]
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

    def _find_data_source(self, data_source_short_name ):
        data_source = self.atomic_db.session.query(DataSource).filter(DataSource.short_name == data_source_short_name).first()
        return

    def _find_atom(self, atomic_number):
        return  self.atomic_db.session.query(Atom).filter(Atom.atomic_number == atomic_number).first()

    def _find_ion(self, atomic_number, ion_number, data_source_short_name= None):
        atom = self._find_atom(atomic_number)
        if data_source_short_name is None:
            return  self.atomic_db.session.query(Ion).filter(Ion.ion_number == ion_number, Ion.atom == atom ).first()
        else:
            data_source = self._find_data_source(data_source_short_name)
            return Ion.query.filter_by(ion_number = ion_number, atom = atom, data_source=data_source).first()


