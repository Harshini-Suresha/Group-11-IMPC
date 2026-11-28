CREATE DATABASE IF NOT EXISTS DCDM;
USE DCDM;

-- DROP TABLE IF EXISTS Parameter;
CREATE TABLE Parameter(
	parameter_id varchar(100) NOT NULL,
	parameter_name varchar(150) NOT NULL,
	parameter_group varchar(50) NOT NULL,
	CONSTRAINT param_pk PRIMARY KEY (parameter_id)
)




DROP TABLE IF EXISTS IMPC_Param;
CREATE TABLE IMPC_Param(
	impcParameterOrigId int NOT  NULL,
    parameterID varchar(100) NOT NULL,
    name varchar(250) NOT NULL,
    description varchar(450) NULL,
    
    CONSTRAINT IMPC_Param_pk PRIMARY KEY (impcParameterOrigId),
    CONSTRAINT param_impc_param_fk FOREIGN KEY (parameterID) REFERENCES Parameter (parameter_id)
)


    

-- DROP TABLE IF EXISTS IMPC_Proc;
CREATE TABLE IMPC_Proc(
	impcParameterOrigId int NOT NULL,
	name varchar(250),
	description varchar(600),
	isMandatory BOOLEAN NOT NULL,
	CONSTRAINT IMPC_Proc_pk PRIMARY KEY (impcParameterOrigId)
)


-- DROP TABLE IF EXISTS Gene;
CREATE TABLE Gene(
	gene_accession_id varchar(15) NOT NULL,
	gene_symbol varchar(30) NOT NULL,
	CONSTRAINT gene_pk PRIMARY KEY (gene_accession_id)
)




-- DROP TABLE IF EXISTS Analysis;
CREATE TABLE Analysis(
	analysis_id varchar(30),
	gene_accession_id varchar(15) NOT NULL,
	mouse_strain varchar(10),
	mouse_life_stage varchar(30),
	parameter_id varchar(100) NOT NULL,
	pvalue double NOT NULL,
	CONSTRAINT analysis_pk PRIMARY KEY (analysis_id),
	CONSTRAINT gene_fk FOREIGN KEY (gene_accession_id) REFERENCES Gene(gene_accession_id),
	CONSTRAINT param__analysis_fk FOREIGN KEY (parameter_id) REFERENCES `Parameter` (parameter_id)
)


-- Display rows that contain 'image' in the 'parameter_group' column
SELECT * FROM `Parameter` WHERE parameter_group like '%image%';



-- Create a new grouping called 'Image'
UPDATE `Parameter`
SET parameter_group = 'Image'
WHERE parameter_name LIKE '%image%'
  AND parameter_group = 'other';
