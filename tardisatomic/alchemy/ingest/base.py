from tardisatomic.alchemy import Ion
class IngestError(ValueError):
    pass


class BaseIngest(object):
    def __init__(self, atomic_db):
        self.atomic_db = atomic_db
        if not self.requirements_satisfied():
            raise IngestError('Requirements for ingest not satisfied')


    def existing_ions(self):
        if self.atomic_db.session.query(Ion).count() > 0:
            return 0


