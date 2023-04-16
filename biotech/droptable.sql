/* drop  table  if exists temp1;
create table temp1 as
select s.SourceIDas SourceID, Min(m1.EValue)  as Evalue
    from matches  as m1
    inner join sequences  as s on s.Accession= m1.Accession
    group by s.SourceID;
/* drop table if exists temp2; 
create table temp2 as
select s.Accessionas Accession, s.SourceIDas SourceID, m2.Evalue  as Evalue
    from matches  as m2
    inner join sequences  as s on s.Accession=m2.Accession;
/* drop table if exists temp3;
create table temp3 as
select t2.Accession as Accession, t1.SourceID as SourceID, t1.Evalue  as Evalue
    from temp1  as t1
    inner join temp2  as t2  on t1.SourceID = t2.SourceID and  t1.Evalue=t2.EValue;
select Name,  Accession, Evalue
from temp3  as t3
inner join organisms as o on t3.SourceID=o.SourceID;
