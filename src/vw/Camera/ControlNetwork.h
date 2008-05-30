

ControlMeasure
Used to record a coordinate (measurement) in an image that corresponds to a control point

- set line/sample
- set line/sample sigma
- set camera model
- set date/time
- set choosername
- set ignore on/off


ControlPoint
- set cp id
- add/remove/access ControlMeasure
- size (returns number of measures)
- set ignore on/off
- type (ground point, or tie point)  <-- or two subclasses?
- set lon/lat/radius
- lon/lat/radius sigma


ControlNet
- set target body name
- set network name/description
- contol point access/size()
- camera model access (per control point)
- add/delete points
- read/write controlnet file
- various methods for computing error
- assoc. with image list/serial number


Things needed for the LM matrices:

- for a given point and measure, return the imaged pixel location
- for a given mesaure, we need the sigma
- position and sigma for gcp
- sigma for position/pose (sigma should come from cmodel?)
