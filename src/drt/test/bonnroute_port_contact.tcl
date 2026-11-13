# BonnRoute DR keeps a FlexDR pin contact on every port. Each port (met3 at
# the left die edge, met2 at the bottom edge, met4 at the top edge; most of
# them half a pitch off track) has the tracks around it blocked next to the
# pin on its layer and the layers below, so the route leaves the pin with
# the access point's own via up, one wire on the layer above, and a via
# down (a bridge). Collapsed into a wrong-way wire on the pin layer, the
# bridge would start inside the pin with an extended end and no connected
# via; FlexDR's connectivity check would then stop with DRT-0206
# (ibex_sky130hd with BonnRoute GR, port data_wdata_o[14]).
# bonnroute_port_contact.{def,guide} come from opus-kit/portfix
# gen_port_contact.py.
source "helpers.tcl"
read_lef "sky130hd/sky130hd.tlef"
read_lef "sky130hd/sky130hd_std_cell.lef"
read_def "bonnroute_port_contact.def"
read_guides "bonnroute_port_contact.guide"

set_bonnroute_profile 2013-improved
detailed_route -bonnroute -verbose 1

puts pass
