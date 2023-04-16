select Name, Accession, Evalue
from
    (select i2.Accession, i1.SourceID, i1.Evalue
    from 
        (select s.SourceID, Min(m1.EValue) as Evalue
        from matches as m1 
        inner join sequences as s on s.Accession= m1.Accession 
        group by s.SourceID) as i1
    inner join 
        (select s.Accession, s.SourceID, m2.Evalue
        from matches as m2 
        inner join sequences as s on s.Accession=m2.Accession) as i2
    on i1.SourceID = i2.SourceID and i1.Evalue=i2.EValue) as i3
inner join organisms as o
on i3.SourceID=o.SourceID
