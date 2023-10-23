library(deSolve)
#easiest to make grid sized as close to the diameter of the roo. big as the the diameter of the root:

# the root is attached to the x-pane of the glass
# y is horizontal depth
# z is vertical depth
setwd("C:\\Users/sgerber/OneDrive - University of Florida/Documents/Students/Amanda/RI-model-test/root-injection-test/")
source("tendency.R")


#geometry and grid
#########################

r0 = 2.5 #enter radius of root

L  =  195.5#length of box
W  =  35
H  =  346 #vertical depth

#temporal grid setup

time = seq(0,3600,1200)

#spatial grid setup in a way that the box size is the smallest possible size to 
# the root cross section, and that the root is completely in one box

#the root is along the z axis, in the middle of the x-axis, at the edge of the y-axis


nx = L %/% (2*r0) # make sure that a single cell is larger than the root dimameter
nx = ifelse((nx %% 2) == 0, nx-1,nx) #make sure it is odd, so that the root in the middle is completely in 1 box
dx = L/nx
x  = seq(0.5*dx,L-0.5*dx,dx)

#same for y direction:
ny = W %/% (2*r0)
dy = W/ny
y  = seq(0.5*dy, W-0.5*dy,dy)

#z direction
nz = round(H/(2*r0))# no constraints needed, just make sure grid size is approximately the sam
dz = H/nz
z = seq(0.5*dz, H-0.5*dz,dz)

#soil
porewater = 0.5*0.7 # [unitless] fraction of volume where diffusion can occur
BD = 1.1 #bulk density ug mm-3

### root setup ######
#####################

Lroot = 95  # root length change number [mm]
H0    =  80 #upper depth at which root is placed [mm]
H1    = Lroot + H0 #lower depth

#calculate where doc is added
weight = numeric(nz)  #fraction of root length within z grid box
mid    =  (nx+1)/2   #index in x direction where root is -> in the middle
top    = seq(0,nz-1,1)*dz #indicate the top and bottom location of z direction
bottom = seq(1,nz,1)*dz   #    ""                   ""

#
#cells that are filled top to bottom with root
idz = which((top>=H0) & (bottom<=H1),arr.ind=TRUE) 
weight[idz] = 1

idzH0 = idz[1]-1   #cell at the top with root part
weight[idzH0] = (bottom[idzH0]-H0)/dz #the fraction of root length within grid

idzH1 = idz[length(idz)] + 1 # cell at the bottom with root part
weight[idzH1]  = (H1 - top[idzH1])/dz #raction of root length within grid

V_straw = r0^2*pi*weight*dz    #volume of root/straw in each box, just a vector
V_free  = array(data=porewater*dx*dy*dz,dim=c(nx,ny,nz)) #free volume to diffuse -> pore space
V_free[mid,1,] = (dx*dy*dz - V_straw)*porewater #special volume in cells where the root is


#estimate parameters for biology
################################
secperyear= 86400*365
Km_amp = 1  # relationship between eq. biomass and KM Km = Km_amp*M
basemin = 0.02/15/secperyear #base mineralization [ug mm-3] -+> #2percent carbon, 15years turnover

#base mineralization 1.1 ug CO2 hr-1 kg-1 -> Shelby
basemin = 1.1/1000 * BD *12/44 /3600  # Shelby's base respiration [CO2 hr-1 kg-1]  #
#d13 -20  

death = 0.033/secperyear  #deathrate [s-1]
cue  = 0.5  #carbon use efficiency [unitless]

## Derive a few parameter to get into the ballpark using equilibrium
# of pre-addition soil. 

#doc uptake coefficient estimated so that M=C at equilibrium
uptake = death^2*(1-cue)/cue^2/basemin #[s-1 (ug mm^3)^-1]

#depolymerization rate so that it is baserelease under equilibrium
Vprime = basemin*(1+Km_amp)/cue  #[ug mm^-3 s-1]   

#Km calculated based on equilibrium biomass
Km = Km_amp*cue*basemin/death/(1-cue) #[ug mm^-3 s-1]

Diffusion_coef = 1e-6 #[mm^2 s-1]

#parameter list not recommended to modify here, but rather change above
pars=c(D=Diffusion_coef,
       porewater=porewater,
       cue=cue,uptake=uptake,
       death=death,
       basemin=basemin) 


#initial conditions
c0 = array(data=0,dim=c(nx,ny,nz)) #doc concentration [ugC mm-3]
m0 = c0 #microbial concentration [ugC mm-3]

#equilibrium glucose concentration [ugC mm-3]
c0 = pars["death"]/pars["uptake"]/pars["cue"]

#input from root exudates
input = 40.1 # amount of exudates added [ug]
flux_in = weight*input*dz/L  #g in each cell
c0[mid,1,] = c0[mid,1,0] + flux_in/V_free[mid,1,] #concentration bump in each cell [ug mm-3]
#c0[mid,1,] = 1 #right now just placeholder


m0[] = pars["basemin"]*pars["cue"]/pars["death"] #equilibrium with base mineralization
times=seq(0,3600,600)

#initial condition vector -> note that ode.3d requires concatenated vectors,
#rather than arrays
state0 = c(as.vector(c0),as.vector(m0))

#just a test to see if tendencies can be calculated
a = dcdt(0,state0,dx,dy,dz,pars,V_free)
stop()

result=ode.3D(y=state0,func=dcdt,times=times,nspec=2,dimens=c(nx,ny,nz),method='rk4',dx=dx,dy=dy,dz=dz,parms=pars,V_free=V_free)

