from astropy import units as u

from sqlalchemy.ext.declarative import declarative_base

from sqlalchemy import Column

from sqlalchemy.orm import relationship

from sqlalchemy import Integer, Float, String, Boolean, ForeignKey

from sqlalchemy import UniqueConstraint


Base = declarative_base()



class Atom(Base):
    __tablename__ = "atoms"
    id = Column(Integer, primary_key=True)
    atomic_number = Column(Integer)
    name = Column(String(150))
    symbol = Column(String(5))
    group = Column(Integer)
    period = Column(Integer)

    def __repr__(self):
        return "{0} Z={1}".format(self.name, self.atomic_number)

class Isotope(Base):
    __tablename__ = "isotopes"
    id = Column(Integer, primary_key=True)
    atom_id = Column(Integer, ForeignKey("atoms.id"))
    mass_number_A = Column(Integer)
    stable = Column(Boolean)


    atom = relationship("Atom", uselist=False, backref='isotopes')


class Decay(Base):
    __tablename__ = "decays"
    id = Column(Integer, primary_key=True)
    source_isotope_id = Column(Integer, ForeignKey("isotopes.id"))
    target_isotope_id = Column(Integer, ForeignKey("isotopes.id"))
    decay_type_id = Column(Integer, ForeignKey("decay_types.id"))
    probability = Column(Float)
    half_life = Column(Float)

    decay_type = relationship("DecayType")

class DecayType(Base):
    __tablename__ = "decay_types"
    id = Column(Integer, primary_key=True)
    name = Column(String(150))
    comment = Column(String(150))




class Ion(Base):
    __tablename__ = "ions"
    id = Column(Integer, primary_key=True)
    atom_id = Column(Integer, ForeignKey("atoms.id"))
    data_source_id = Column(Integer, ForeignKey("data_sources.id"))
    ion_number = Column(Integer)
    ionization_energy = Column(Float)
    ionization_energy_unit = u.eV

    data_source = relationship('DataSource', backref='ions')
    atom = relationship("Atom", uselist=False, backref='ions')

    def __repr__(self):
        return "{0} Ion={1}".format(self.__tablename__, self.ion_number)


class Level(Base):
    __tablename__ = "levels"

    id = Column(Integer, primary_key=True)
    ion_id = Column(Integer, ForeignKey("ions.id"))
    level_number = Column(Integer)
    g = Column(Integer)
    energy = Column(Float)
    label = Column(String(120))
    theoretical = Column(Boolean)
    data_source_id = Column(Integer, ForeignKey('data_sources.id'))


    ion = relationship("Ion", uselist=False, backref='levels')
    data_source = relationship('DataSource', backref='levels')

    #__table_args__ = (UniqueConstraint('level_number', 'label' ,'energy', 'ion_id' ,'data_source_id',  name='_level_uc'),
    #                 )

    def __repr__(self):
        return "{0} level_number={1}".format(self.__tablename__, self.level_number)

class TransitionType(Base):
    __tablename__ = "transition_types"
    id = Column(Integer, primary_key=True)
    name = Column(String(150))
    comment = Column(String(150))

class TransitionValueType(Base):
    __tablename__ = "transition_value_types"
    id = Column(Integer, primary_key=True)
    name = Column(String(150))
    unit_id = Column(Integer, ForeignKey('units.id'))
    dtype_id = Column(Integer, ForeignKey("dtypes.id"))

    dtype = relationship('DType')
    unit = relationship('Unit')


class TransitionValue(Base):
    __tablename__ = "transition_values"
    id = Column(Integer, primary_key=True)
    transition_id = Column(Integer, ForeignKey("transitions.id"))
    transition_value_type_id = Column(Integer, ForeignKey(
        "transition_value_types.id"))
    value = Column(String(150))
    comment = Column(String(150))

    data_source_id = Column(Integer, ForeignKey("data_sources.id"))

    transition_value_type = relationship('TransitionValueType', backref='transition_values')
    data_source = relationship('DataSource', backref='transition_values')
    transition = relationship('Transition', backref='transition_value')


class Transition(Base):
    __tablename__ = "transitions"
    id = Column(Integer, primary_key=True)
    transition_type_id = Column(Integer, ForeignKey("transition_types.id"))
    source_level_id = Column(Integer, ForeignKey("levels.id"))
    target_level_id = Column(Integer, ForeignKey("levels.id"))

    data_source_id = Column(Integer, ForeignKey("data_sources.id"))
    comment = Column(String(150))

    transition_type = relationship("TransitionType", uselist=False,
                                   backref='transitions')

    source_level = relationship(
        "Level", primaryjoin=(Level.id==source_level_id), uselist=False)
    target_level = relationship(
        "Level", primaryjoin=(Level.id==target_level_id), uselist=False)

    data_source = relationship('DataSource', backref='transitions')


class DataSource(Base):
    __tablename__ = "data_sources"
    id = Column(Integer, primary_key=True)
    short_name = Column(String(20))
    name = Column(String(120))
    description = Column(String(800))
    data_source_quality = Column(Integer)



class DType(Base):
    __tablename__ = "dtypes"
    id = Column(Integer, primary_key=True)
    short_name = Column(String(10))
    name = Column(String(150))

class Unit(Base):
    __tablename__ = 'units'

    id = Column(Integer, primary_key=True)
    unit = Column(String(150))

