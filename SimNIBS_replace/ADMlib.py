#!/usr/local/bin/python
import os
import logging
import numpy as np
import fmm3dpy as fmm
import time
import numpy.linalg as la

logger = logging.getLogger('simnibs')
#*********************************************************************
#  This library compute E-fields induced during TMS via reciprocity and auxiliary
#  dipoles. The method is implemented for magnetic dipole sources.
#  companion pre-print: https://www.biorxiv.org/content/10.1101/2020.05.27.120022v2
#  created by Luis Gomez (July 2020)!
#  Last modified by Xu Zhang (July 2024)
#*********************************************************************

'''
   Authors of auxiliary dipole method (ADM) and code for determining rapidly the
   optimum coil position and orientation:
   Luis J. Gomez, Moritz Dannhauer, and Angel V. Peterchev;
   Duke University, Durham, North Carolina, U.S.A.

   The development of the Duke ADM algorithm and code were supported by the
   National Institute of Mental Health and the National Institute on Aging of the
   National Institutes of Health under Award Numbers K99MH120046, RF1MH114268,
   RF1MH114253, and U01AG050618. The content is solely the responsibility of the authors
   and does not necessarily represent the official views of the National Institutes of Health.
   
   The copyrights of this software are owned by Duke University. As such, two licenses to this software are offered:

      * An open source license under the GNU General Public License (GPL-v2.0) (https://opensource.org/licenses/gpl-2.0.php).
      * A custom license with Duke University, for use without the GPL-v2.0 restrictions. 

   As a recipient of this software, you may choose which license to receive the code under.
   Outside contributions to the Duke owned code base cannot be accepted unless the contributor transfers
   the copyright to those changes over to Duke University.

   To enter a license agreement without the GPL-v2.0 restrictions,
   please contact the Digital Innovations department at Duke Office of Licensing and Ventures
   (https://olv.duke.edu/software/) at olvquestions@duke.edu with reference to “OLV File No. 7148” in your email. 

   Please note that this software is distributed AS IS, WITHOUT ANY WARRANTY;
   and without the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
   
   Modified by Xu Zhang, 2024 (each change is indicated in the script)
    
'''

def recipcode(rv,jv,rs,ks,A):
	#this function computes E-fields via reciprocity no auxiliary dipoles
	#rv is 3 by ntetra and has mesh tetrahedron centroid positions
	#jv is 3 by ntetra and has total conduction current at each tetrahedron
	#rs is 3 by ncoil and has the coil dipole positions centered about the origin
	#ks is 3 by ncoil and has the coil dipole weights
	#A is 4 by 4 by number of coilpositions and each 4 by 4 matrix is a translation of the coil to a point above the scalp

	#paramameter:
	prec=10**(-3); #this number determines the accuracy of H-primary evaluation higher accuracy=slower time increases ~log(accuracy)

	#reads in sizes of each array
	npos=A.shape[2]; #number of scalp positions
	ncoil=rs.shape[1]; 
	nsource=rv.shape[1];
	ntarget=ncoil*npos;
	
	#generate copies of coil
	rp=np.transpose(np.c_[np.transpose(rs) ,np.ones([ncoil])]); #pads coil positions to 4 by ncoil
	robs=np.array(np.zeros([3,ntarget]));
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		robs[:,st:en]=np.matmul(A[0:3,:,i],rp);

	start = time.time()
	logger.debug("Computing H-primary");
	Hprimary=computeHprimary(rv,jv,robs,prec);
	end = time.time()
	logger.debug(f"H-primary time: {end-start:.2f}s");
	kp=ks;
	Etotal=np.zeros(npos);
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		kp=np.matmul(A[0:3,0:3,i],ks);
		Etotal[i]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
	return Etotal;

# Modified by Xu Zhang, 2024 (start)
def ADM(rv,jv,rs,ks,A,coildir,resampleflag,N,Xr,Yr,Zr):
    #this function computes E-fields via reciprocity with auxiliary dipoles
    #rv is 3 by ntetra and has mesh tetrahedron centroid positions
    #jv is 3 by ntetra and has total conduction current at each tetrahedron
    #rs is 3 by ncoil and has the coil dipole positions centered about the origin
    #ks is 3 by ncoil and has the coil dipole weights
    #A is 4 by 4 by number of coilpositions and each 4 by 4 matrix is a translation of the coil to a point above the scalp
    #coildir is a 3 by number of coil orientations that gives the y-direction orientation
    #N is a 1-by-3 array that indicates the number of Gauss points along each axis

    #paramameter:
    prec=10**(-3); #this number determines the accuracy of H-primary evaluation higher accuracy=slower time increases ~log(accuracy)
    # N=[17,17,2]; % Modified by Xu Zhang, 2024
    #generate auxiliary dipoles
    Nj=coildir.shape[1];
    
    # Modified by Xu Zhang, 2024 (start)
    if resampleflag:
        raux,kaux=resamplecoil(rs,ks,N,Nj,coildir,Xr,Yr,Zr);
    else:
        raux = rs;
        kaux = ks[:,:,np.newaxis];
    # Modified by Xu Zhang, 2024 (end)
    
    #reads in sizes of each array
    npos=A.shape[2]; #number of scalp positions
    ncoil=raux.shape[1]; 
    nsource=rv.shape[1];
    ntarget=ncoil*npos;
    
    #generate copies of coil
    rp=np.transpose(np.c_[np.transpose(raux) ,np.ones([ncoil])]); #pads coil positions to 4 by ncoil
    robs=np.array(np.zeros([3,ntarget]));
    for i in range(0,npos):
        st=i*ncoil;
        en=(i+1)*ncoil;
        robs[:,st:en]=np.matmul(A[0:3,:,i],rp);

    start = time.time()
    logger.debug("Computing H-primary");
    Hprimary=computeHprimary(rv,jv,robs,prec);
    end = time.time()
    logger.debug(f"H-primary time: {end-start:.2f}s")
    kp=kaux[:,:,0];
    
    Etotal=np.zeros([Nj,npos]);
    for i in range(0,npos):
        st=i*ncoil;
        en=(i+1)*ncoil;
        for j in range(0,Nj):
            kp=np.matmul(A[0:3,0:3,i],kaux[:,:,j]);#for flat coils you only need the 3rd element
            Etotal[j,i]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
    return Etotal;
# Modified by Xu Zhang, 2024 (end)

def recipcodemag(rv,jvx,jvy,jvz,rs,ks,A):
	#this function computes E-field unidirectional approximation of the magnitude via reciprocity no auxiliary dipoles
	#rv is 3 by ntetra and has mesh tetrahedron centroid positions
	#jv is 3 by ntetra and has total conduction current at each tetrahedron
	#rs is 3 by ncoil and has the coil dipole positions centered about the origin
	#ks is 3 by ncoil and has the coil dipole weights
	#A is 4 by 4 by number of coilpositions and each 4 by 4 matrix is a translation of the coil to a point above the scalp

	#paramameter:
	prec=10**(-3); #this number determines the accuracy of H-primary evaluation higher accuracy=slower time increases ~log(accuracy)

	#reads in sizes of each array
	npos=A.shape[2]; #number of scalp positions
	ncoil=rs.shape[1]; 
	nsource=rv.shape[1];
	ntarget=ncoil*npos;
	
	#generate copies of coil
	rp=np.transpose(np.c_[np.transpose(rs) ,np.ones([ncoil])]); #pads coil positions to 4 by ncoil
	robs=np.array(np.zeros([3,ntarget]));
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		robs[:,st:en]=np.matmul(A[0:3,:,i],rp);
	Etotal=np.zeros([npos,3]);
	kp=ks;
	start = time.time()
	logger.debug("Computing H-primary");
	Hprimary=computeHprimary(rv,jvx,robs,prec);
	end = time.time()	
	logger.info(f"H-primary time: {end-start:.2f}s")
	start = end
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		kp=np.matmul(A[0:3,0:3,i],ks);
		Etotal[i,0]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
	start = time.time()
	logger.debug("Computing H-primary");
	Hprimary=computeHprimary(rv,jvy,robs,prec);
	end = time.time()	
	logger.info(f"H-primary time: {end-start:.2f}s")
	start = end
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		kp=np.matmul(A[0:3,0:3,i],ks);
		Etotal[i,1]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
	start = time.time()
	logger.debug("Computing H-primary");
	Hprimary=computeHprimary(rv,jvz,robs,prec);
	end = time.time()	
	logger.info(f"H-primary time: {end-start:.2f}s")
	start = end
	for i in range(0,npos):
		st=i*ncoil;
		en=(i+1)*ncoil;
		kp=np.matmul(A[0:3,0:3,i],ks);
		Etotal[i,2]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
	Etotal=np.sqrt(Etotal[:,0]**2+Etotal[:,1]**2+Etotal[:,2]**2);
	return Etotal;

# Modified by Xu Zhang, 2024 (start)
def ADMmag(rv,jvx,jvy,jvz,rs,ks,A,coildir,resampleflag,N,Xr,Yr,Zr):
    #this function computes E-field unidirectional approximation of the magnitude via reciprocity with auxiliary dipoles
    #rv is 3 by ntetra and has mesh tetrahedron centroid positions
    #jv is 3 by ntetra and has total conduction current at each tetrahedron
    #rs is 3 by ncoil and has the coil dipole positions centered about the origin
    #ks is 3 by ncoil and has the coil dipole weights
    #A is 4 by 4 by number of coilpositions and each 4 by 4 matrix is a translation of the coil to a point above the scalp
    #coildir is a 3 by number of coil orientations that gives the y-direction orientation
    #N is a 1-by-3 array that indicates the number of Gauss points along each axis

    #paramameter:
    prec=10**(-3); #this number determines the accuracy of H-primary evaluation higher accuracy=slower time increases ~log(accuracy)
    # N=[17,17,2]; # Modified by Xu Zhang, 2024
    #generate auxiliary dipoles
    Nj=coildir.shape[1];
    
    # Modified by Xu Zhang, 2024 (start)
    if resampleflag:
        raux,kaux=resamplecoil(rs,ks,N,Nj,coildir,Xr,Yr,Zr);
    else:
        raux = rs;
        kaux = ks[:,:,np.newaxis];
    # Modified by Xu Zhang, 2024 (end)
    
    #reads in sizes of each array
    npos=A.shape[2]; #number of scalp positions
    ncoil=raux.shape[1]; 
    nsource=rv.shape[1];
    ntarget=ncoil*npos;
    
    #generate copies of coil
    rp=np.transpose(np.c_[np.transpose(raux) ,np.ones([ncoil])]); #pads coil positions to 4 by ncoil
    robs=np.array(np.zeros([3,ntarget]));
    for i in range(0,npos):
        st=i*ncoil;
        en=(i+1)*ncoil;
        robs[:,st:en]=np.matmul(A[0:3,:,i],rp);

    start = time.time()
    Etotal=np.zeros([Nj,npos,3]);
    logger.debug("Computing H-primary");
    Hprimary=computeHprimary(rv,jvx,robs,prec);
    end = time.time()
    logger.info(f"H-primary time: {end-start:.2f}s")
    start = end
    kp=kaux[:,:,0];
    for i in range(0,npos):
        st=i*ncoil;
        en=(i+1)*ncoil;
        for j in range(0,Nj):
            kp=np.matmul(A[0:3,0:3,i],kaux[:,:,j]);#for flat coils you only need the 3rd element
            Etotal[j,i,0]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
    Hprimary=computeHprimary(rv,jvy,robs,prec);
    end = time.time()
    logger.info(f"H-primary time: {end-start:.2f}s")
    start = end
    for i in range(0,npos):
        st=i*ncoil;
        en=(i+1)*ncoil;
        for j in range(0,Nj):
            kp=np.matmul(A[0:3,0:3,i],kaux[:,:,j]);#for flat coils you only need the 3rd element
            Etotal[j,i,1]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
    Hprimary=computeHprimary(rv,jvz,robs,prec);
    end = time.time()
    logger.info(f"H-primary time: {end-start:.2f}s")
    start = end
    for i in range(0,npos):
        st=i*ncoil;
        en=(i+1)*ncoil;
        for j in range(0,Nj):
            kp=np.matmul(A[0:3,0:3,i],kaux[:,:,j]);#for flat coils you only need the 3rd element
            Etotal[j,i,2]=-np.inner(Hprimary[:,st:en].flatten(),kp.flatten());
    Etotal=np.sqrt(Etotal[:,:,0]**2+Etotal[:,:,1]**2+Etotal[:,:,2]**2);
    return Etotal;
# Modified by Xu Zhang, 2024 (end)

def computeHprimary(rs,js,robs,prec):
	#this function computes H-fields via FMM3D library
	#convention is 3 by number of points
	#prec determines the accuracy of the multipole expansion
	#Note: for magnetic dipoles electromagnetic duality implies that
	#if we pass magnetic dipoles weights as js we get negative E-primary.
	# As such, this function is used to compute E-primary due to magnetic currents also.
	outex=fmm.Output();
	Hprimary=np.zeros([3,robs.shape[1]]);
	muover4pi=-1e-7;
	js1=js[0,:];
	out=fmm.lfmm3d(eps=prec,sources=rs,targets=robs,charges=js1,pgt=2);
	logger.debug("Run 1")
	Hprimary[1,:]=-out.gradtarg[2,:];
	Hprimary[2,:]=out.gradtarg[1,:];
	js1=js[1,:];
	out=fmm.lfmm3d(eps=prec,sources=rs,targets=robs,charges=js1,pgt=2);
	logger.debug("Run 2")
	Hprimary[0,:]=out.gradtarg[2,:];
	Hprimary[2,:]=Hprimary[2,:]-out.gradtarg[0,:];
	js1=js[2,:];
	out=fmm.lfmm3d(eps=prec,sources=rs,targets=robs,charges=js1,pgt=2);
	logger.debug("Run 3")
	Hprimary[0,:]=Hprimary[0,:]-out.gradtarg[1,:];
	Hprimary[1,:]=Hprimary[1,:]+out.gradtarg[0,:];
	Hprimary *= muover4pi
	return Hprimary;

# Modified by Xu Zhang, 2024 (start)
def resamplecoil(rs,ks,N,Nj,coildir,Xr,Yr,Zr):
    #rs is 3 by ncoil and has the coil dipole positions centered about the origin
    #ks is 3 by ncoil and has the coil dipole weights
    #N is 3 by 1 and has the number of auxiliary dipoles along each dimension
    #Nj is an integer number of orientations 
    #coildir is 3 by Nj and has the y orientation of the coil
    #create copies of coil with different orientations
    rs2=np.zeros([3,rs.shape[1],Nj]);
    ks2=np.zeros([3,ks.shape[1],Nj]);
    #x=y cross z = y[1] x -y[0] y
    for i in range(0,Nj):
        rs2[2,:,i]=rs[2,:];
        rs2[2,:,i]=rs[2,:];
        rs2[0,:,i]= rs[0,:]*coildir[1,i]+rs[1,:]*coildir[0,i];
        rs2[1,:,i]=-rs[0,:]*coildir[0,i]+rs[1,:]*coildir[1,i];
        ks2[2,:,i]=ks[2,:];
        ks2[0,:,i]= ks[0,:]*coildir[1,i]+ks[1,:]*coildir[0,i];#unnecessary ks is z oriented
        ks2[1,:,i]=-ks[0,:]*coildir[0,i]+ks[1,:]*coildir[1,i];#unnecessary ks is z oriented
    #find range for interpolation
    
    # Modified by Xu Zhang, 2024 (start)
    if len(Xr) == 0:
        Xm=min(rs2[0,:,:].flatten());
        Xp=max(rs2[0,:,:].flatten());
    else:
        Xm = Xr[0];
        Xp = Xr[1];
    
    if len(Yr) == 0:
        Ym=min(rs2[1,:,:].flatten());
        Yp=max(rs2[1,:,:].flatten());
    else:
        Ym = Yr[0];
        Yp = Yr[1];
    
    if len(Zr) == 0:
        Zm=min(rs2[2,:,:].flatten());
        Zp=max(rs2[2,:,:].flatten());
            #this assumes the original coil has 3 layers and the quadrature is gaussian
            #will need to be changed for other coils
        Zd=.6127016653792583*0.5*(Zp-Zm);
        Zm=Zm-Zd;
        Zp=Zp+Zd;
    else:
        Zm = Zr[0];
        Zp = Zr[1];
    # Modified by Xu Zhang, 2024 (end)
    
    #get gauss quadrature nodes for interpolation from txt file upto N=50
    interpnodes = np.loadtxt(os.path.join(os.path.dirname(__file__), 'gauss50.txt'))
    #translate quadrature rules to coil box
    XX=0.5*(interpnodes[int(N[0]*(N[0]-1)/2):int(N[0]*(N[0]+1)/2)]+1.0)*(Xp-Xm)+Xm;
    YY=0.5*(interpnodes[int(N[1]*(N[1]-1)/2):int(N[1]*(N[1]+1)/2)]+1.0)*(Yp-Ym)+Ym;
    ZZ=0.5*(interpnodes[int(N[2]*(N[2]-1)/2):int(N[2]*(N[2]+1)/2)]+1.0)*(Zp-Zm)+Zm;
    #generate auxiliary grid
    Y,Z,X=np.meshgrid(YY,ZZ,XX);
    raux=np.zeros([3,X.size]);
    kaux=np.zeros([3,X.size,Nj]);
    raux[0,:]=X.flatten();
    raux[1,:]=Y.flatten();
    raux[2,:]=Z.flatten();
    del X;
    del Y;
    del Z;
    #generate auxiliary dipole weights
    for kk in range(0,Nj):
        Lx=lagrange(rs2[0,:,kk],XX);
        Ly=lagrange(rs2[1,:,kk],YY);
        Lz=lagrange(rs2[2,:,kk],ZZ);
        for k in range(0,N[2]):
            for j in range(0,N[1]):
                for i in range(0,N[0]):
                    L=Lx[i,:]*Ly[j,:]*Lz[k,:];
                    kaux[0,int(i+(j+k*N[1])*N[0]),kk]=np.inner(ks2[0,:,kk],L);
                    kaux[1,int(i+(j+k*N[1])*N[0]),kk]=np.inner(ks2[1,:,kk],L);
                    kaux[2,int(i+(j+k*N[1])*N[0]),kk]=np.inner(ks2[2,:,kk],L);
    return raux,kaux
# Modified by Xu Zhang, 2024 (end)

def lagrange(x,pointx):
	n=pointx.size;
	L=np.ones([n,x.size]);
	for i in range(0,n):
		for j in range(0,n):
			if i != j:
				L[i,:]=np.multiply(L[i,:],(x-pointx[j])/(pointx[i]-pointx[j]));
	return L;

if __name__ == '__main__':
    test_computeHprimary()
