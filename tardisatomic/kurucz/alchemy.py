from sqlalchemy.ext.declarative import declarative_base

from sqlalchemy import Column

from sqlalchemy import Integer, Float, String

Base = declarative_base()

class Lines(Base):
    id = Column(Integer, primary_key=True)
