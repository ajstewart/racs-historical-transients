# Generated by Django 2.2.5 on 2020-02-19 10:46

from django.db import migrations


class Migration(migrations.Migration):

    initial = False

    dependencies = [
        ('images', '0001_initial')
    ]

    operations = [
        migrations.RunSQL(
            ["CREATE EXTENSION IF NOT EXISTS q3c;"],#upgrade
            ["DROP EXTENSION IF EXISTS q3c;"],#downgrade
        ),
        migrations.RunSQL(
            ["CREATE INDEX ON images_crossmatches (q3c_ang2ipix(ra, dec));"],
            ["DROP INDEX images_crossmatches_q3c_ang2ipix_idx;"],
        ),
        migrations.RunSQL(
            ["CLUSTER images_crossmatches_q3c_ang2ipix_idx ON images_crossmatches;"],
            [],
        ),
        migrations.RunSQL(
            ["ANALYZE images_crossmatches;"],
            [],
        ),
    ]
