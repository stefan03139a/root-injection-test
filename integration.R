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

L  =  195.5#length of box [mm]
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
basemin = 0.02/15/secperyear #base mineralization [ug mm-3] -+> #2percent carbon, 15years turnover


#basic values to derive parameters 
#base mineralization 1.1 mg CO2 hr-1 kg-1 -> Shelby
#converting into ugC mm-3 s-1
basemin = 1.1e-3 * BD *12/44 /3600  
Km_amp = 1  # relationship between eq. biomass and KM Km = Km_amp*M
lambda_c = 1/86400  #turnover of DOC is assumed 1 day-1

#direct parameters
death = 12/secperyear  #deathrate [s-1]
cue  = 0.5  #carbon use efficiency [unitless]
f_doc = 0.3  #fraction of microbial death become doc
doubling_time = 3600*4  #maximum microbial growth rate expressed in time for 2x

#derived parameters and state variables equilibria
umax = log(2)/(cue*doubling_time)
Vprime = basemin*(1+Km_amp) #maximum depolymerization rate
m0 = basemin*cue/death/(1-f_doc*cue) #equilibrum microbes
c0 = death/cue/uptake
Km = Km_amp*m0  #microbial half saturation for depoly 
Kc = umax/lambda_c*m0 -c0 #half saturation constant for doc

#water
#https://bionumbers.hms.harvard.edu/bionumber.aspx?id=104089&ver=7
#600 um^2/s
#0.6e-5 cm2/s = 0.6e-3 mm2/s  = 0.6e3 um2/s
D=0.6e-3  #conversion into mm2/s-1

#Chenu and Roberson 0.5e-6 cm2/s @ 65 % WFPS
D=0.5e-4  #conversion int mm2/s

D=6.3e-4*0.28*0.7 #Priesack and Kisser-Priesack, 1993, f=0.28 (impedance at 100% water)

#parameter list not recommended to modify here, but rather change above
pars=c(D=Diffusion_coef,
       porewater=porewater,
       cue=cue,umax=umax,
       f_doc,
       death=death,
       Km = Km,
       Kc=Kc,
       Vprime=Vprime) 


#initial conditions
c0 = array(data=c0,dim=c(nx,ny,nz)) #doc concentration [ugC mm-3]
m0 = array(data=m0,dim=c(nx,ny,nz))

#input from root exudates
input = 1.68e6  # amount of exudates added [ug]
flux_in = weight*input*dz/L  #ug in each cell

c0[mid,1,] = c0[mid,1,]+ flux_in/V_free[mid,1,]

times=seq(0,3600,600)

#initial condition vector -> note that ode.3d requires concatenated vectors,
#rather than arrays
state0 = c(as.vector(c0),as.vector(m0))

#just a test to see if tendencies can be calculated
#a = dcdt(0,state0,dx,dy,dz,pars,V_free)
#stop()

result=ode.3D(y=state0,func=dcdt,times=times,nspec=2,dimens=c(nx,ny,nz),method='rk4',dx=dx,dy=dy,dz=dz,parms=pars,V_free=V_free)



#########################
#Evaluations
#########################


#test array
data=matrix(data=1:81,nrow=3)
dim(result)
nt=length(time)
doc = array(data=result[,2:(nx*ny*nz+1)],dim=c(nt,nx,ny,nz))
m   = array(data=result[,(nx*ny*nz+2):(2*nx*ny*nz+1)],dim=c(nt,nx,ny,nz))

plot(time,m[,nx/2,ny/2,mid])
