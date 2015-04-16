from sqlalchemy.ext.declarative import declarative_base

from sqlalchemy import Column

from sqlalchemy.orm import relationship

from sqlalchemy import Integer, Float, String, Boolean, ForeignKey

Base = declarative_base()



class Atoms(Base):
    __tablename__ = "Atoms"
    id = Column(Integer, primary_key=True)
    atomic_number = Column(Integer)
    name = Column(String(150))
    atomic_sym = Column(String(5))
    group = Column(Integer)
    period = Column(Integer)

    ions = relationship("Ion")
    isotope = relationship("Isotope")


class Isotope(Base):
    __tablename__ = "Isotope"
    id = Column(Integer, primary_key=True)
    atom_id = Column(Integer, ForeignKey("Atoms.id"))
    mass_number_A = Column(Integer)
    stable = Column(Boolean)
    decay = relationship("Decay")


class Decay(Base):
    __tablename__ = "Decay"
    id = Column(Integer, primary_key=True)
    source_isotope_id = Column(Integer, ForeignKey("Isotope.id"))
    target_isotope_id = Column(Integer, ForeignKey("Isotope.id"))
    decay_type_id = Column(Integer, ForeignKey("Decay_type.id"))
    probability = Column(Float)
    half_life = Column(Float)


class DecayTyp(Base):
    __tablename__ = "Decay_type"
    id = Column(Integer, primary_key=True)
    name = Column(String(150))
    comment = Column(String(150))

    decay = relationship("Decay")


class Ion(Base):
    __tablename__ = "Ion"
    id = Column(Integer, primary_key=True)
    atom_id = Column(Integer, ForeignKey("Atoms.id"))
    ion_number = Column(Integer)
    ionization_energy_at_ground_level = Column(Float)

    transition = relationship("Transition")
    level = relationship("Level")


class Level(Base):
    __tablename__ = "Level"

    id = Column(Integer, primary_key=True)
    ion_id = Column(Integer, ForeignKey("Ion.id"))
    level_number = Column(Integer)
    g = Column(Integer)
    energy = Column(Float)
    ionization_energy = Column(Float)
    label = Column(String(120))
    theoretical = Column(Boolean)
    data_source_id = Column(Integer, ForeignKey('Data_source.id'))

    transition = relationship("Transition")


class TransitionType(Base):
    __tablename__ = "Transition_type"
    id = Column(Integer, primary_key=True)
    name = Column(String(150))
    comment = Column(String(150))

    transitions = relationship("Transition")


class TransitionValueType(Base):
    __tablename__ = "Transition_value_type"
    id = Column(Integer, primary_key=True)
    name = Column(String(150))
    unit = Column(String(150))
    dtype_id = Column(Integer, ForeignKey("Dtype.id"))

    transition_type = relationship("Transition_value")


class TramsitionValue(Base):
    __tablename__ = "Transition_value"
    id = Column(Integer, primary_key=True)
    transition_id = Column(Integer, ForeignKey("Transition.id"))
    transition_value_type_id = Column(Integer, ForeignKey("Transition_value_type.id"))
    value = Column(String(150))
    comment = Column(String(150))

    data_source_id = Column(Integer, ForeignKey("Data_source.id"))


class Transition(Base):
    __tablename__ = "Transition"
    id = Column(Integer, primary_key=True)
    transition_type_id = Column(Integer, ForeignKey("Transition_type.id"))
    source_level_id = Column(Integer, ForeignKey("Level.id"))
    target_level_id = Column(Integer, ForeignKey("Level.id"))
    source_ion_id = Column(Integer, ForeignKey("Ion.id"))
    target_ion_id = Column(Integer, ForeignKey("Ion.id"))

    data_source_id = Column(Integer, ForeignKey("Data_source.id"))
    comment = Column(String(150))

    transition_value = relationship("Transition_value")


class DataSource(Base):
    __tablename__ = "Data_source"
    id = Column(Integer, primary_key=True)
    short_name = Column(String(20))
    name = Column(String(120))
    description = Column(String(800))
    data_source_quality = Column(Integer)

    transition_value = relationship("Transition_value")
    transition = relationship("Transition")
    level = relationship("Level")


class DType(Base):
    __tablename__ = "Dtype"
    id = Column(Integer, primary_key=True)
    short_name = Column(String(10))
    name = Column(String(150))

    TransitionValueType = relationship("Transition_value_type")


