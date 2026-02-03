## Sweep Configuration 

Sets up parameter sweeps for the electrostatics solver 

Also implements Nelder Mead optimization of metrics computed after simulation



TODOs:

Mesh based CIV use: FindPointsGSLIB

It is supposed to be quicker than the way I am working at the moment - but need to look into it


### Optimization


#### CIV

This entire code block is very messy due to prolongued debugging. 

The CIV is computed form the cathode wires top edge upwards. Two modes are available random sampling and informed sweep, where in the first n random samples are drawn (subject to sampling error - also not really feasible due to computational complexity) and in the latter we try to isolate regions of charge insensitive volume with the knowledge that if a volume is charge insensitive everything to its right will be. This is provided the electron terminates on the wall, given the height definition due to field leakage some electrons will travel downwards and Hit the cathode plane, in this case we apply the same logic downwards. 

###### Random sample

This works, but it produces few truly charge insensitive seeds, typically we capture only seeds in areas of high mesh element density due to their prominance. 


###### informed sweep

This does not work - it terminates after finding one or two points as it cant find any CI neighbors

TODO CLean up code and try again

Here we utilize a binary search starting at max_r - geom_tol and half way up the tpc we then half the volume and check again, this is repeated until the first charge insensitive point is found, this is done batched to save on time so we evaluate the half quarter eight etc points in parallel. Once it is found the boundary is identified, from this point on we move down and towards lower r, repeating the procedure row wise until all charge insensitve mesh elements at the boundary to charge sensitive are found. Lastly the row at the bottom will be isolated. Once all boundaries are known a precomputed index of mesh element center coordinates is used to mark all elements below the bottom plane (where traces hit the cathode) and elements to the bottom right of the elements idenfied as hitting the wall are marked as charge insensitive and finally the mesh centric volume is computed.  

###### Column sweep

The only reliable version at the moment

Still requires that we do some form of gridding, basically this does bisection on columns until a charge insensitive element is found and then it proceeds to find the highest. These are then used to form a boundary and all volumes below this are summed as charge insensitive
