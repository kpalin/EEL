DROP TABLE IF EXISTS sequence;
DROP TABLE IF EXISTS region;
DROP TABLE IF EXISTS cismodule;
DROP TABLE IF EXISTS sites;
DROP TABLE IF EXISTS tfactor;
DROP TABLE IF EXISTS pairs;
DROP TABLE IF EXISTS alnCols;
DROP TABLE IF EXISTS gfffile;





drop table if exists sequence;
drop table if exists region;
drop table if exists cismodule;
drop table if exists tfactor;
drop table if exists sliceCisLink;
drop table if exists slice;
DROP TABLE IF EXISTS alnCols ;

drop table if exists sites;
drop table if exists gene_description;

create table gene_description(
	enscode varchar(40),
	modified date,
	description text,
        KEY (enscode),
	FOREIGN KEY (enscode) REFERENCES slice(enscode)
);

CREATE TABLE sequence (
	id INT AUTO_INCREMENT PRIMARY KEY,
	species VARCHAR(125),
	DBname VARCHAR(125) NOT NULL,
	name VARCHAR(125) NOT NULL,
	UNIQUE (name,species,DBname)
);


-- mysqlimport -p  bbu_palin -v --local gff.txt

CREATE TABLE slice (
  id int(11) AUTO_INCREMENT PRIMARY KEY,
  beginPos int(11) NOT NULL default '0',
  endPos int(11) NOT NULL default '0',
  fileName varchar(200) NOT NULL default '',
  enscode varchar(40),
	geneName varchar(20),
  seqName varchar(200),
  revComplement int(11) NOT NULL default 0,
	seqID INT,
	KEY (seqID),
  KEY `geneIndex` (`enscode`),
  KEY (fileName,seqName)
);


create table cismodule (
	id INT AUTO_INCREMENT PRIMARY KEY,
	score FLOAT,
	tmpID varchar(255),
	KEY (tmpID)
);

	
CREATE TABLE sliceCisLink (
	sliceID int not null,
	cisID int not null,
	KEY (sliceID,cisID),
	PRIMARY KEY (cisID,sliceID)
);


CREATE TABLE alnCols (
	id INT AUTO_INCREMENT PRIMARY KEY,
	tfID INT,
	scoredelta FLOAT,
	cisID INT NOT NULL,
	cisColId INT NOT NULL,
	KEY (cisID,cisColId),
	KEY (tfID)
);

CREATE TABLE sites (
	id INT AUTO_INCREMENT PRIMARY KEY,
	weight FLOAT,
	pos INT,
	strand CHAR(1),
	regID INT NOT NULL,
	colID INT NOT NULL,
	KEY (regID),
	KEY (colID)
);



CREATE TABLE region (
	id INT AUTO_INCREMENT PRIMARY KEY,
	seqID INT,
	beginPos INT,
	endPos INT,
	cisID INT,
	KEY (seqID),
	KEY (cisID)
);


CREATE TABLE tfactor (
	id INT AUTO_INCREMENT PRIMARY KEY,
	width INT,
	name VARCHAR(255) UNIQUE,
	KEY (name)
);




