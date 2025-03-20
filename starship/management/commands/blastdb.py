from django.core.management.base import BaseCommand, CommandError
from tempfile import TemporaryDirectory
from starship import models as starship_models
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
from starship import views
import subprocess
from django.conf import settings


class Command(BaseCommand):
    help = 'creates blast database for protein and nucleotide'

    def handle(self, *args, **options):
        self.phage_blastdb()
        self.annotation_blastdb()

    # creates nucleotide blast database
    def phage_blastdb(self):
        # temp = TemporaryDirectory()
        # file_path = os.path.join(temp.name, "AMD_Complete_Phage_Genomes")
        file_path = settings.NUCLEOTIDE_DATABASE + '.fa'
        phages = starship_models.Starship.objects.all()
        phage_list = []
        for phage in phages:
            sequence = SeqRecord(Seq(starship.starship_sequence, generic_nucleotide), id=starship.starship_name,
                                 description=starship.starship_name)
            phage_list.append(sequence)

        SeqIO.write(phage_list, file_path, "fasta")
        # nucleotide_db = os.path.join(settings.NUCLEOTIDE_DATABASE, 'AMD_Complete_Phages')
        subprocess.run(['makeblastdb', '-in', '"%s"' % file_path, '-input_type', 'fasta', '-dbtype', 'nucl', '-title',
                        'MAS Starships', '-out', settings.NUCLEOTIDE_DATABASE], check=True)

    # creates sequence blast database
    def annotation_blastdb(self):
        # temp = TemporaryDirectory()
        # file_path = os.path.join(temp.name, "AMD_Annotated_Proteins.fasta")
        file_path = settings.PROTEIN_DATABASE + '.fa'
        annotations = starship_models.Annotation.objects.all()
        annotation_list = []
        for annotation in annotations:
            anno = annotation.annotation
            aa_sequence = annotation.sequence
            public_note = annotation.public_notes
            private_note = annotation.private_notes
            flag = annotation.get_flag_display()
            phages = views.annotation_starships(annotation)
            sequence = SeqRecord(Seq(aa_sequence, generic_protein), id=annotation.accession + " |",
                                 description="%s | %s | %s | %s %s" % (anno, public_note, private_note, flag, phages))
            annotation_list.append(sequence)
        SeqIO.write(annotation_list, file_path, "fasta")
        # protein_db = os.path.join(settings.PROTEIN_DATABASE, 'AMD_Annotated_Proteins.faa')
        subprocess.run([
            'makeblastdb',
            '-in', '"%s"' % file_path,
            '-input_type', 'fasta',
            '-dbtype', 'prot',
            '-parse_seqids',
            '-title', 'MAS Proteins',
            '-out', settings.PROTEIN_DATABASE
        ], check=True)
