#!/bin/sh

#
# $Log$
#

echo Fetching descriptions from the Sanger center (ENSEMBL).

mysql -h kaka.sanger.ac.uk -u anonymous  homo_sapiens_core_19_34b -B -e "select stable_id,modified,description  from gene_stable_id s, gene_description d  where s.gene_id=d.gene_id and modified<'2004-01-01'" >tmp$$.tab


echo Loading data to local database (gene_description table).

mysql -u bbu-palin -ppalindrome bbu_palin<<EOF
drop table if exists gene_description;

create table gene_description(
	enscode varchar(40),
	modified date,
	description text,
        KEY (enscode),
	FOREIGN KEY (enscode) REFERENCES slice(enscode)
);

LOAD DATA LOCAL INFILE 'tmp$$.tab' INTO TABLE gene_description;

DELETE FROM gene_description WHERE enscode='stable_id';
EOF

echo Clean up.
rm tmp$$.tab
echo Done.
