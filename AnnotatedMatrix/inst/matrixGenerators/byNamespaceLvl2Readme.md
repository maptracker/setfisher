### Template file

This is a template file that is copied into the second-level
subdirectories of the `byNamespace` directories. This text (above the
horizontal rule) will not be copied.

----
# Namespace Pair <NS1> + <NS2>

This folder holds matrix files with both the `<NS1>` and `<NS2>`
namespaces, irrespective of whether one is in rows or columns (that
is, <NS1> could be associated with either the rows or columns of any
given matrix here). The files are actually symbolic links back to the
[byAuthority][BA] directory.

### CAUTION CAUTION CAUTION

Files may be loaded directly from this folder tree, but their location
is not persistent! The tree's primary use is to aid in exploration of
available files. If you are encoding a file in a production workflow,
please reference instead the `byAuthority` file paths represented by
the targets of each symlink. Those locations are static and will
persist as new versions are loaded.

[BA]: ../../../byAuthority
