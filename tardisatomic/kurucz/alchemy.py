from sqlalchemy.ext.declarative import declarative_base

from sqlalchemy import Column

from sqlalchemy.orm import relationship, backref

from sqlalchemy import Integer, Float, String, Boolean, ForeignKey, Table

Base = declarative_base()

level_line_association_table = Table('level_line_association', Base.metadata,
    Column('Levels_id', Integer, ForeignKey('Levels.id')),
    Column('Lines_id', Integer, ForeignKey('Lines.id'))
)


class Atoms(Base):
    __tablename__ = "Atoms"
    id = Column(Integer, primary_key=True)
    atomic_number = Column(Integer)
    name = Column(String(150))
    atomic_sym = Column(String(5))
    group = Column(Integer)
    period = Column(Integer)

    ions = relationship("Ions")




class Ions(Base):
    __tablename__ = "Ions"
    id = Column(Integer, primary_key=True)
    atom_id = Column(Integer, ForeignKey("Atoms.id"))
    ion_number = Column(Integer)
    ionization_energy_at_ground_level = Column(Float)

    levels = relationship("Levels")




class Levels(Base):
    __tablename__ = "Levels"

    id = Column(Integer, primary_key=True)
    ion_id = Column(Integer, ForeignKey("Ions.id"))
    level_number = Column(Integer)
    g = Column(Integer)
    energy = Column(Float)
    ionization_energy = Column(Float)
    label = Column(String(120))
    theoretical = Column(Boolean)
    data_source_id = Column(Integer, ForeignKey('Data_Source.id'))
    dataSource = relationship('Data_Source')


class Lines(Base):
    __tablename__ = "Lines"
    id = Column(Integer)
    level_id_lower = relationship("Levels", secondary=level_line_association_table, backref="Lines")
    level_id_upper = relationship("Levels", secondary=level_line_association_table, backref="Lines")
    wavelength = Column(Float)
    f_ul = Column(Float)
    f_lu = Column(Float)
    loggf = Column(Float)


class DataSource(Base):
    __tablename__ = "Data_Source"
    id = Column(Integer, primary_key=True)
    short_name = Column(String(20))
    name = Column(String(120))
    description = Column(String(800))
    data_source_quality = Column(Integer)



