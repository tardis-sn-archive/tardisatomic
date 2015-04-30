from tardisatomic.alchemy import Ion, Atom, Level,Transition, TransitionType,\
    TransitionValue ,TransitionValueType, DataSource

from sqlalchemy.orm import  aliased
import pandas as pd


class CreateHd5(object):
    def __init__(self, atomic_db, hd5_filename):
        self.atomic_db = atomic_db


    def load_data(self):
        pass


    def _load_atoms(self):
        atomic_data = self.atomic_db.session.query(Atom).order_by(
            Atom.atomic_number).all()
        self.atom_data = self._to_data_frame(atomic_data)

    def _load_ions(self):
        ion_data = self.atomic_db.session.query(Ion, Atom).join("atom").values(
            Atom.atomic_number, Ion.ion_number,
            Ion.ionization_energy)
        self.ion_data = self._to_data_frame(ion_data)

    def _load_transitions(self):
        target_level = aliased(Level)
        source_level = aliased(Level)
        trans_a = aliased(Transition)
        trans_b = aliased(Transition)
        trans_data = self.atomic_db.session.query(TransitionValue,
                                                  TransitionValue.value,
                                                  Transition,
                                                  target_level.id.label(
                                                      "target_level_id"),
                                                  source_level.id.label(
                                                      "source_level_id"),
                                                  target_level.level_number.label(
                                                      "target_level_number"),
                                                  source_level.level_number.label("source_level_number"),
                                                  TransitionType.name,
                                                  TransitionValueType.name
                                                  ).join(Transition).join(
            target_level,
            Transition.target_level_id == target_level.id).join(source_level,
            Transition.target_level_id == source_level.id).join(
            TransitionType).filter(TransitionType.name == 'Line')
        bla = self._to_data_frame(trans_data)
        a = 1
        pass
        #transition_data =


    def _load_levels(self):
        level_data = self.atomic_db.session.query(Level, Ion, Atom).join(
            "ion", "atom"). values(Atom.atomic_number, Ion.ion_number,
                                   Level.level_number, Level.g, Level.energy)
        self.level_data = self._to_data_frame(level_data)

    def _to_data_frame(self, sql_data):
        rec_data = [rec.__dict__ for rec in sql_data]
        return pd.DataFrame.from_records(rec_data)