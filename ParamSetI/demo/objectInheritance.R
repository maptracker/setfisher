## Example of how ParamSetI can be inheritted by your own
## ReferenceClass objects

## A toy object class is defined in another demo - I hear you like
## demos, so I put demos in your demos:
demo("exampleParamSetObject", package="ParamSetI", ask=FALSE)

## Create a new object, set a custom value for the "inc" parameter
pse <- ParamSetExample( params=list(inc=10L) )

## Poke it a few times
pse$inc()
pse$dec()
pse$inc()

## Display it:
pse

## Check all parameters
pse$showParameters()

## Change a parameter:
pse$param("increment", 25L)
pse$inc()
pse



