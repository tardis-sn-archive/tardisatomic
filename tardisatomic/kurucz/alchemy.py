from sqlalchemy.ext.declarative import declarative_base

from sqlalchemy import Column

from sqlalchemy.orm import relationship, backref

from sqlalchemy import Integer, Float, String, Boolean, ForeignKey

Base = declarative_base()

class Atoms(Base):
    id = Column(Integer, primary_key=True)
    atomic_number = Column(Integer)


class Ions(Base):
    __tablename__ = "Ions"
    id = Column(Integer, primary_key=True)
    atomic_number = Column(Integer)
    ion_number = Column(Integer)

    levels = relationship("Levels")




class Levels(Base):
    __tablename__ = "Levels"

    id = Column(Integer, primary_key=True)

    atomic_number = Column(Integer)
    ion_number = Column(Integer)
    level_number = Column(Integer)
    g = Column(Integer)
    energy = Column(Float)
    label = Column(String(120))
    theoretical = Column(Boolean)
    data_source_id = Column(Integer, ForeignKey('Data_Source.id'))
    dataSource = relationship('Data_Source')

class Lines(Base):
    __tablename__ = "Lines"
    id = Column(Integer)
    atomic_number = Column(Integer)
    ion_number = Column(Integer)
    level_number = Column(Integer)


class DataSource(Base):
    __tablename__ = "Data_Source"
    id = Column(Integer, primary_key=True)
    short_name = Column(String(20))
    name = Column(String(120))
    description = Column(String(800))
    data_source_quality = Column(Integer)



