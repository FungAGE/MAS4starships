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
        self.stdout.write('Creating BLAST databases...')
        self.starship_blastdb()
        self.annotation_blastdb()
        self.swissprot_blastdb()

    # creates nucleotide blast database
    def starship_blastdb(self):
        self.stdout.write('Creating Starship nucleotide BLAST database...')
        file_path = settings.NUCLEOTIDE_FASTA_PATH
        starships = starship_models.Starship.objects.all()
        starship_list = []
        for starship in starships:
            sequence = SeqRecord(Seq(starship.starship_sequence), id=starship.starship_name,
                               description=starship.starship_name)
            starship_list.append(sequence)

        SeqIO.write(starship_list, file_path, "fasta")
        subprocess.run(['makeblastdb', '-in', '"%s"' % file_path, '-input_type', 'fasta', '-dbtype', 'nucl', '-title',
                      'MAS Starships', '-out', settings.NUCLEOTIDE_DB_PATH], check=True)
        self.stdout.write(self.style.SUCCESS('Successfully created Starship BLAST database'))

    # creates sequence blast database
    def annotation_blastdb(self):
        self.stdout.write('Creating internal protein BLAST database...')
        file_path = settings.PROTEIN_FASTA_PATH
        annotations = starship_models.Annotation.objects.all()
        annotation_list = []
        for annotation in annotations:
            anno = annotation.annotation
            aa_sequence = annotation.sequence
            public_note = annotation.public_notes
            private_note = annotation.private_notes
            flag = annotation.get_flag_display()
            starships = views.annotation_starships(annotation)
            sequence = SeqRecord(Seq(aa_sequence), id=annotation.accession + " |",
                               description="%s | %s | %s | %s %s" % (anno, public_note, private_note, flag, starships))
            annotation_list.append(sequence)
        SeqIO.write(annotation_list, file_path, "fasta")
        subprocess.run([
            'makeblastdb',
            '-in', '"%s"' % file_path,
            '-input_type', 'fasta',
            '-dbtype', 'prot',
            '-parse_seqids',
            '-title', 'MAS Proteins',
            '-out', settings.PROTEIN_DB_PATH
        ], check=True)
        self.stdout.write(self.style.SUCCESS('Successfully created internal protein BLAST database'))

    # creates swissprot blast database
    def swissprot_blastdb(self):
        self.stdout.write('Creating SwissProt BLAST database...')
        swissprot_fasta = settings.SWISSPROT_FASTA_PATH
        swissprot_db = settings.SWISSPROT_DB_PATH
        
        # Check if source file exists
        if not os.path.exists(swissprot_fasta):
            self.stdout.write(self.style.WARNING(
                f'SwissProt FASTA file not found at {swissprot_fasta}. Skipping SwissProt database creation.'
            ))
            return

        # Create the BLAST database
        try:
            subprocess.run([
                'makeblastdb',
                '-in', swissprot_fasta,
                '-dbtype', 'prot',
                '-parse_seqids',
                '-title', 'SwissProt',
                '-out', swissprot_db
            ], check=True)
            self.stdout.write(self.style.SUCCESS('Successfully created SwissProt BLAST database'))
        except subprocess.CalledProcessError as e:
            self.stdout.write(self.style.ERROR(f'Failed to create SwissProt BLAST database: {str(e)}'))
