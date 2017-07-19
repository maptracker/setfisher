## ParamSetI - A ReferenceClass interface for parameter management

This is an R ReferenceClass object (code in [./ParamSetI][source],
tar.gz in [../packageReleases][targz]) designed to be inherited
(`contains`) by other RefClass objects. Parameters are stored as
key/value pairs, where the keys are simple strings and the values are
generally (but not exclusively) single 'scalars' (generally a string
or number).

It handles get/set operations, specification of definitions for
parameter keys, and reporting of current values and definitions to the
user.

Imports [CatMisc][CatMisc] and [EventLogger][EventLogger].

[source]: ParamSetI
[targz]: ../packageReleases

[CatMisc]: https://github.com/maptracker/CatMisc/
[EventLogger]: https://github.com/maptracker/setfisher/EventLogger
