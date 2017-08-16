### Template file

This is a template file that is copied into the `byAuthority`
directories. This text (above the horizontal rule) will not be copied.

----
# AnnotatedMatrix Files - Primary folder

This folder holds [MatrixMarket][MTX] files that can be read and
manipulated by the [AnnotatedMatrix][AM] R package. It is organized by
the data authority - the primary source of the data used to generate
the files. Each folder on this level is an authority, and each
subfolder is a release version. Within each subfolder will be one or
more MatrixMarket files. For example:

```
byAuthority/
  HGNC/
    2017-07-26/
      Map@EntrezGene_to_EnsemblGene@HGNC@2017-07-26.mtx
      Map@HGNC_to_EnsemblGene@HGNC@2017-07-26.mtx
      Map@HGNC_to_EntrezGene@HGNC@2017-07-26.mtx
      Map@HGNC_to_IMGT@HGNC@2017-07-26.mtx
      ... etc ...
    2017-10-30/
      Map@EntrezGene_to_EnsemblGene@HGNC@2017-10-30.mtx
      Map@HGNC_to_EnsemblGene@HGNC@2017-10-30.mtx
      Map@HGNC_to_EntrezGene@HGNC@2017-10-30.mtx
      Map@HGNC_to_IMGT@HGNC@2017-10-30.mtx
      ... etc ...
    
```

These files are the 'primary' data files. The `byAuthority` folder
will have sibling folders that organize the files in different
directory hierarchies. Those hierarchies will not contain files, but
rather symlinks pointing back to the 'primary' files in this
folder. These symlink hierarchies are provided to aid human browsing
of the matrix files.

In addition to the `.mtx` flat files, the directory may also contain
`.mtx.rds` files. These are serialized R data structures created after
the `.mtx` file is first parsed. In general, you should continue to
work directly with the `.mtx` files - the code can detect updates
(based on date stamps) and refresh the RDS files as needed.


[MTX]: http://math.nist.gov/MatrixMarket/formats.html
[AM]: https://github.com/maptracker/setfisher/tree/master/AnnotatedMatrix
