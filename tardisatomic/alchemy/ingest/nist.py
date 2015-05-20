from tardisatomic.alchemy.ingest import BaseIngest
from tardisatomic.nist.io import download_ionization, parse_ionization_data
from tardisatomic.alchemy import Ion, Atom

class NISTIonization(BaseIngest):
    requirements = []

    def requirements_satisfied(self):
        return True


    def ingest(self, spectra='h-uuo'):
        nist_download = download_ionization(spectra)
        ionization_data = parse_ionization_data(nist_download)

        for i, row in ionization_data.iterrows():
            ion = Ion(ion_number=row.ion_number,
                      ionization_energy=row.ionization_energy)
            ion.atom = self.atomic_db.session.query(Atom).filter_by(
                atomic_number=row.atomic_number).one()
            self.atomic_db.session.add(ion)

        self.atomic_db.session.commit()





