CREATE DATABASE IF NOT EXISTS gr11;
USE gr11;


CREATE TABLE mouse (
mouse_id int AUTO_INCREMENT PRIMARY KEY,
strain varchar(5) NOT NULL, 
CONSTRAINT min_length_strain CHECK (CHAR_LENGTH(strain)>=3),
life_stage varchar(17) NOT NULL,
CONSTRAINT min_length_life_stage CHECK (CHAR_LENGTH(life_stage)>=4)
);

CREATE TABLE mouse_mgi (
mouse_mgi_id varchar(11) PRIMARY KEY NOT NULL,
CONSTRAINT min_length_mouse_mgi_id CHECK (CHAR_LENGTH(mouse_mgi_id)>=9)
);

CREATE TABLE original_parameter (
original_parameter_id int PRIMARY KEY,
original_description varchar(500),
original_name varchar(255)
);


CREATE TABLE parameter (
parameter_id varchar(20) PRIMARY KEY,
CONSTRAINT min_length_parameter_id CHECK (CHAR_LENGTH(parameter_id)>=15),
parameter_description varchar(500),
parameter_name varchar(74) NOT NULL,
CONSTRAINT min_length_parameter_name CHECK (CHAR_LENGTH(parameter_name)>=2),
original_parameter_id int UNIQUE,
    CONSTRAINT fk_original_parameter
    FOREIGN KEY (original_parameter_id)
    REFERENCES original_parameter(original_parameter_id)
);


CREATE TABLE analysis (
analysis_id varchar(15) PRIMARY KEY,
CONSTRAINT min_length_analysis_id CHECK (CHAR_LENGTH(analysis_id)>=15),
p_value double NOT NULL,
CONSTRAINT check_p_value CHECK (p_value BETWEEN 0 and 1),
mouse_mgi_id varchar(11) NOT NULL,
    CONSTRAINT fk_mouse_mgi1
    FOREIGN KEY (mouse_mgi_id)
    REFERENCES mouse_mgi(mouse_mgi_id),
parameter_id varchar(20),
    CONSTRAINT fk_parameter1
    FOREIGN KEY (parameter_id)
    REFERENCES parameter(parameter_id),
mouse_id int,
	CONSTRAINT fk_mouse
	FOREIGN KEY (mouse_id)
	REFERENCES mouse(mouse_id)
);




CREATE TABLE omim (
omim_id varchar(11) PRIMARY KEY NOT NULL
);


CREATE TABLE do_disease (
do_disease_id varchar(12) PRIMARY KEY NOT NULL,
disease_name varchar(255)
);


CREATE TABLE mouse_mgi_gene_symbol (
mouse_mgi_gene_symbol_id int AUTO_INCREMENT PRIMARY KEY,
gene_symbol varchar(13) NOT NULL,
CONSTRAINT min_length_gene_symbol CHECK (CHAR_LENGTH(gene_symbol)>=1),
mouse_mgi_id varchar(11) NOT NULL,
CONSTRAINT fk_mouse_mgi2
    FOREIGN KEY (mouse_mgi_id)
    REFERENCES mouse_mgi(mouse_mgi_id)
);


CREATE TABLE mouse_mgi_omim_do_disease (
mouse_mgi_omim_do_disease_id int AUTO_INCREMENT PRIMARY KEY,
mouse_mgi_id varchar(11) NOT NULL,
    CONSTRAINT fk_mouse_mgi3
    FOREIGN KEY (mouse_mgi_id)
    REFERENCES mouse_mgi(mouse_mgi_id),
omim_id varchar(11) NOT NULL,
    CONSTRAINT fk_omim
    FOREIGN KEY (omim_id)
    REFERENCES omim(omim_id),
do_disease_id varchar(12) NOT NULL,
	CONSTRAINT fk_do_disease
	FOREIGN KEY (do_disease_id)
	REFERENCES do_disease(do_disease_id)

);

CREATE TABLE parameter_groupings (
parameter_id varchar(20),
    CONSTRAINT fk_parameter2
    FOREIGN KEY (parameter_id)
    REFERENCES parameter(parameter_id),
grouping_name varchar(255), 
grouping_name_id int AUTO_INCREMENT PRIMARY KEY
);


