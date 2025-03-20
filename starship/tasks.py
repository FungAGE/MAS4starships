from tempfile import TemporaryDirectory
import os
import subprocess
from configparser import ConfigParser

from celery import shared_task

# from celery.app.control import Inspect
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from django.db import transaction
from django.conf import settings
from django.db.models.signals import post_save

from starship.genomic_loci_conversions import *
from starship import gene_calling
from starship import models as starship_models
from starship.forms import parse_prots_from_coords

from MAS.celery import app


def on_starship_uploaded_or_removed(
    sender, **kwargs
):  # TODO: ADD SIGNAL WHICH CALLS THIS UPON GENOME DELETION
    if settings.INTERNAL_NUCLEOTIDE_DB_PATH:
        i = app.control.inspect(settings.CELERY_WORKERS)
        queue = i.scheduled()
        if queue:
            for worker in queue.keys():
                for task in queue[worker]:
                    if (
                        task["request"]["name"]
                        == "starship.tasks.create_internal_nucleotide_blastdb"
                    ):
                        app.control.revoke(task["request"]["id"])

            create_internal_nucleotide_blastdb.apply_async(countdown=5 * 60)


def on_annotation_changed(sender, **kwargs):
    i = app.control.inspect(settings.CELERY_WORKERS)
    queue = i.scheduled()
    if queue:
        for worker in queue.keys():
            for task in queue[worker]:
                if (
                    task["request"]["name"]
                    == "starship.tasks.create_internal_protein_blastdb"
                ):
                    app.control.revoke(task["request"]["id"])

        create_internal_protein_blastdb.apply_async(countdown=5 * 60)


starship_models.starship_upload_complete.connect(on_starship_uploaded_or_removed)
post_save.connect(on_annotation_changed, sender=starship_models.Annotation)
starship_models.starship_upload_complete.connect(on_annotation_changed)


@shared_task
def create_internal_nucleotide_blastdb():
    print("Creating internal nucleotide BLAST database.")

    if settings.INTERNAL_NUCLEOTIDE_DB_PATH:
        fasta_file_path = settings.INTERNAL_NUCLEOTIDE_DB_PATH + ".fa"
        starships = starship_models.Starship.objects.all()
        starship_list = []
        for starship in starships:
            sequence = SeqRecord(
                Seq(starship.starship_sequence),
                id=starship.starship_name,
                description=starship.starship_name,
            )
            starship_list.append(sequence)

        SeqIO.write(starship_list, fasta_file_path, "fasta")
        subprocess.run(
            [
                "makeblastdb",
                "-in",
                fasta_file_path,
                "-input_type",
                "fasta",
                "-dbtype",
                "nucl",
                "-title",
                "MAS Starships",
                "-out",
                settings.NUCLEOTIDE_DATABASE,
            ],
            check=True,
        )


@shared_task
def create_internal_protein_blastdb():
    print("Creating internal protein BLAST database.")

    # Get path from luigi config
    cfg = ConfigParser()
    cfg.read(settings.LUIGI_CFG)
    db_path = cfg["Blastp"]["internal"]
    fasta_file_path = db_path + ".fa"

    annotations = starship_models.Annotation.objects.all()
    annotations = annotations.prefetch_related("feature_set")
    annotations = annotations.prefetch_related("feature_set__starship")

    annotation_list = []
    for annotation in annotations:
        anno = annotation.annotation
        aa_sequence = annotation.sequence
        public_note = annotation.public_notes
        private_note = annotation.private_notes
        flag = annotation.get_flag_display()
        genomes = ", ".join(
            {feature.starship.starship_name for feature in annotation.feature_set.all()}
        )
        record = SeqRecord(
            Seq(aa_sequence),
            id=annotation.accession + " |",
            description="%s | %s | %s | %s | %s"
            % (anno, public_note, private_note, flag, genomes),
        )
        annotation_list.append(record)
    print("writing fasta")
    SeqIO.write(annotation_list, fasta_file_path, "fasta")
    # protein_db = os.path.join(settings.PROTEIN_DATABASE, 'AMD_Annotated_Proteins.faa')
    subprocess.run(
        [
            "makeblastdb",
            "-in",
            '"%s"' % fasta_file_path,
            "-input_type",
            "fasta",
            "-dbtype",
            "prot",
            "-parse_seqids",
            "-title",
            "MAS Internal Protein Database",
            "-out",
            db_path,
        ],
        check=True,
    )


# Create new annotation objects and associated features for CDS
def create_CDS_annotations(
    glimmer_results_file_path, genome, assign_to, new_annotations, new_features
):
    for cds in gene_calling.parse_glimmer_results(glimmer_results_file_path):
        sequence = Seq(starship.starship_sequence)
        protein = get_protein_sequence(cds.start, cds.stop, cds.strand, sequence)

        if starship_models.Annotation.objects.filter(sequence=protein).count() > 0:
            new_annotation = starship_models.Annotation.objects.get(sequence=protein)
        elif protein in new_annotations:
            new_annotation = new_annotations[protein]
        else:
            new_annotation = starship_models.Annotation()
            new_annotation.sequence = protein
            new_annotation.assigned_to = assign_to
            new_annotations[protein] = new_annotation

        new_feature = starship_models.Feature()
        new_feature.genome = genome
        new_feature.start = cds.start
        new_feature.stop = cds.stop
        new_feature.type = "CDS"
        new_feature.strand = cds.strand
        new_feature.annotation = new_annotation
        new_features.append(new_feature)


# Create an annotation and feature objects for given coordinate file
def create_custom_CDS_annotations(
    coordinate_file, translation_table, genome, assign_to, new_annotations, new_features
):
    starship_sequence = Seq(starship.starship_sequence)

    for protein_sequence, cds in parse_prots_from_coords(
        coordinate_file, starship_sequence, translation_table
    ):

        if (
            starship_models.Annotation.objects.filter(sequence=protein_sequence).count()
            > 0
        ):
            new_annotation = starship_models.Annotation.objects.get(
                sequence=protein_sequence
            )

        elif protein_sequence in new_annotations:
            new_annotation = new_annotations[protein_sequence]

        else:
            new_annotation = starship_models.Annotation()
            new_annotation.sequence = protein_sequence
            new_annotation.assigned_to = assign_to
            new_annotations[protein_sequence] = new_annotation

        new_feature = starship_models.Feature()
        new_feature.genome = genome
        new_feature.start = cds.start
        new_feature.stop = cds.stop
        new_feature.type = "CDS"
        new_feature.strand = cds.strand
        new_feature.annotation = new_annotation
        new_features.append(new_feature)


# Create new annotation objects and associated features for tRNAs
def create_trna_annotations(
    trnascan_results_file_path, genome, assign_to, new_annotations, new_features
):
    for tRNA in gene_calling.parse_trnascan_results(trnascan_results_file_path):
        sequence = Seq(starship.starship_sequence)
        rna = get_rna_sequence(tRNA.start, tRNA.stop, tRNA.strand, sequence)

        if starship_models.Annotation.objects.filter(sequence=rna).count() > 0:
            new_annotation = starship_models.Annotation.objects.get(sequence=rna)
        elif rna in new_annotations:
            new_annotation = new_annotations[rna]
        else:
            new_annotation = starship_models.Annotation()
            new_annotation.annotation = "tRNA-%s-%s" % (tRNA.amino, tRNA.codon)
            new_annotation.public_notes = (
                "Called with tRNAscan-SE, Inf score = %f" % tRNA.score
            )
            new_annotation.private_notes = (
                "This annotation was automatically generated."
            )
            new_annotation.sequence = rna
            new_annotation.flag = 8
            new_annotation.assigned_to = assign_to
            new_annotations[rna] = new_annotation

        new_feature = starship_models.Feature()
        new_feature.genome = genome
        new_feature.start = tRNA.start
        new_feature.stop = tRNA.stop
        new_feature.type = "tRNA"
        new_feature.strand = tRNA.strand
        new_feature.annotation = new_annotation
        new_features.append(new_feature)


def add_annotations_and_features_to_db(new_annotations, new_features):
    # Adding objects to db gives them an id
    starship_models.Annotation.objects.bulk_create(list(new_annotations.values()))

    # We must assign newly created id in feature object before creation
    for feat in new_features:
        feat.annotation_id = feat.annotation.id

    starship_models.Feature.objects.bulk_create(new_features)
