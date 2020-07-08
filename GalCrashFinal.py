import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5 import uic
import sys
import pyqtgraph as pg
import time


class Galaxy:
    """A class used to define the initial parameters of a generated galaxy. It also contains functions
    that are used in calculating the evolution of the galaxy. Here we generally deal with the center 
    of the galaxy when looking at parameters such as position."""
    def __init__(self,galmass,ahalo,vhalo,rthalo,galpos,galvel):
        self.galmass = galmass
        self.ahalo = ahalo
        self.vhalo = vhalo
        self.rthalo = rthalo
        self.galpos = galpos
        self.galvel = galvel
        self.galacc = np.full((3,1),0.)

        
    def setPosvel(self,pos,vel):
        self.galpos = pos
        self.galvel = vel
    
    def scaleMass(self,massFact):
        """The Scale Mass Function calculates the parameters of the companion galaxy. This is done by taking the
        ratio of the mass of the companion galaxy to the main galaxy as the input."""
        self.galmass = self.galmass*massFact
        self.vhalo = 1.0*massFact**0.25
        self.ahalo = 0.1*massFact**0.5
        a2 = -self.galmass/(self.vhalo**2)
        a1 = -2.*self.ahalo*self.galmass/(self.vhalo**2)
        a0 = -self.galmass*(self.ahalo**2)/(self.vhalo**2)
        q = a1/3.0 - (a2**2)/9.0
        r = (a1*a2-3.*a0)/6.0 - (a2**3)/27.0

        s1 = (r + np.sqrt((q**3)+(r**2)))**0.333
        s2 = ((r - np.sqrt((q**3)+(r**2)))**0.333)

        self.rthalo=(s1+s2)-a2/3

        


    def MoveGalaxy(self,dtime): 
        """Move Galaxy evolves the position and velocity of the galaxy with input dtime"""
        newpos = self.galpos + self.galvel * dtime + 0.5 * self.galacc *(dtime**2)
        newvel = self.galvel + self.galacc * dtime

        self.galpos = newpos;
        self.galvel = newvel;




    def Acceleration(self,posin):
        """Acceleration function takes the position of the object (star/galaxy) as input and returns the acceleration
         of the object."""
        G = 1.0
        dpos = posin - self.galpos
        #dx = dpos[0]
        #dy = dpos[1]
        #dz = dpos[2]
        r = np.sqrt(np.sum(dpos**2, axis = 0))
        #r = np.sqrt((dx**2)+(dy**2)+(dz**2))
        AccMag= -(G*self.InteriorMass(r))/(r**2)
        calcacc = (dpos*AccMag)/r
        
        return calcacc



    def Potential(self, posin):
        """Potential calculates the potential of the object (star/galaxy) based on the input position given."""
        G=1.0;

        dpos = posin - self.galpos
        #dx = dpos[0]
        #dy = dpos[1]
        #dz = dpos[2]
        r = np.sqrt(np.sum(dpos**2, axis = 0))
        #r = np.sqrt((dx**2)+(dy**2)+(dz**2))
        pot= G*self.InteriorMass(r)/r

        return pot



    def InteriorMass(self, r):
        """Interior Mass returns the mass of the object (star/galaxy) based on its position in relation to the center
        of the galaxy"""
        
        
        indices = r < self.rthalo
        

        intmass = np.full(r.shape, 0.)

        if intmass[indices].shape != (0,):
        #If r< self.rthalo, then we are dealing with a star confined within the radius of the galaxy, and the interior mass
        # is calculated as follows
            intmass[indices] = (self.vhalo**2)*(r[indices]**3)/((self.ahalo+r[indices])**2)
        if intmass[~indices].shape != (0,):
        #Otherwise, the interior mass is simply the mass of the galaxy
            intmass[~indices] = self.galmass
        
        return intmass



    def Density(self, r):
        """Determines the density based on the input r which is the position of the object in relation to the center of
        the galaxy."""
        rinner = r*0.99
        router = r*1.01
        minner = self.InteriorMass(r*0.99)
        mouter = self.InteriorMass(r*1.01)
        dm = (mouter-minner)
        vol=(4/3)*np.pi*((router**3)-(rinner**3))
        density=dm/vol

        return density


    def DynFric(self,pmass,ppos,pvel):
        """Dynamic Friction is used when considering friction in our model. It calculates the friction based on the mass
        position and velocity of the object. The resultant friction contributes to the overall acceleration of that object
        """

        G=1.0
        lnGamma=3.0
        dv = pvel - self.galvel
        v = np.linalg.norm(dv)
        dr = ppos - self.galpos
        r = np.linalg.norm(dr)
        galrho = self.Density(r)
        fricmag = 4.0*np.pi*G*lnGamma*pmass*galrho*v/((1+v)**3)
        friction = (-dv/v)*fricmag
        
        return friction




    def PrintGalaxy(self):
        """Prints information about the properties of the generated galaxy"""
        print("Mass: ", self.galmass)
        print("ahalo:", self.ahalo)
        print("vhalo:", self.vhalo)
        print("rthalo:", self.rthalo)
        print("Position:", self.galpos)
        print("Velocity:", self.galvel)
        print("Acceleration:", self.galacc)



class StarGalaxy(Galaxy):
    """A sub-class of galaxy, used to combine information from the stars and the galactic center to produce a galaxy
    with orbiting stars."""
    def __init__(self,galmass,ahalo,vhalo,rthalo,galpos,galvel,diskSize,galtheta,galphi,n):
        super().__init__(galmass,ahalo,vhalo,rthalo,galpos,galvel)
        self.diskSize = diskSize
        self.galtheta = galtheta*np.pi/180
        self.galphi = galphi*np.pi/180
        self.n = n
        
        #Define the star position, velocity, and acceleration in such a way that they can be called upon later more easily.
        self.starpos = np.full((3,self.n), 0.)
        self.starvel = np.full((3,self.n), 0.)
        self.staracc = np.full((3,self.n), 0.)

    def MoveStars(self,dtime):
        """Function for moving a star over time by calculating the new position and velocity"""
        newstarpos =  self.starpos + self.starvel * dtime + 0.5 * self.staracc * (dtime**2)
        newstarvel = self.starvel + self.staracc * dtime
       
        self.starpos = newstarpos
        self.starvel = newstarvel

    def InitStars(self):
        """InitStars initializes the stars in the galaxy. It generates n number of stars and gives them random positions
        within the disksize of the galaxy. The velocities of the stars are also calculated based on the positions as well as
        the input angles phi and theta"""
        cosphi = np.cos(self.galphi)
        sinphi = np.sin(self.galphi)
        costheta = np.cos(self.galtheta)
        sintheta = np.sin(self.galtheta)
        for i in range(self.n):
            bad = True
            while bad:
                xtry = self.diskSize*(1.-2.*np.random.random())
                ytry = self.diskSize*(1.-2.*np.random.random())
                rtry = np.sqrt(xtry**2+ytry**2)
                if (rtry < self.diskSize): bad = False
            
            ztry = 0.0
            xrot = xtry*cosphi + ytry*sinphi*costheta + ztry*sinphi*sintheta
            yrot = -xtry*sinphi + ytry*cosphi*costheta + ztry*cosphi*sintheta
            zrot = -ytry*sintheta + ztry*costheta
            rot = np.array([xrot,yrot,zrot])
            self.starpos[:,i] = rot + self.galpos.reshape(-1)
            
            vcirc = np.sqrt(self.InteriorMass(rtry)/rtry)

            vxtry = -vcirc*ytry/rtry
            vytry = vcirc*xtry/rtry
            vztry = 0.0

            vxrot = vxtry*cosphi + vytry*sinphi*costheta + vztry*sinphi*sintheta
            vyrot = -vxtry*sinphi + vytry*cosphi*costheta + vztry*cosphi*sintheta
            vzrot = -vytry*sintheta + vztry*costheta

            vrot = np.array([vxrot,vyrot,vzrot])
            self.starvel[:,i] = vrot + self.galvel.reshape(-1)
            self.staracc = np.full((1,3),0.)
    def scaleMass(self, massFact):
        """The function scalemass calculates the disk size of a companion galaxy based on the mass ratio of the
        companion galaxy to the main galaxy"""
        self.diskSize = self.diskSize*np.sqrt(massFact)
        super().scaleMass(massFact)

        
        
class Orbit:
    """The Orbit class calculates initial position and velocity of the two galaxies based on a parabolic orbit"""
    def __init__(self,energy,rp,tp,eccentricity,m1,m2,bod1pos,bod2pos,bod1vel,bod2vel):
        self.energy = energy
        self.rp = rp
        self.tp = tp
        self.eccentricity = eccentricity
        self.m1 = m1
        self.m2 = m2
        self.bod1pos = bod1pos 
        self.bod2pos = bod2pos
        self.bod1vel = bod1vel
        self.bod2vel = bod2vel
        self.initOrbit()

        
        
    def initOrbit(self):
        """InitOrbit initializes the orbit of the two galaxies based on the input masses of the galaxies
        This function assumes a parabolic orbit and returns the position and velocity of both galaxies"""
        #Parabolic Orbit
        mu = self.m1 + self.m2

        p = 2*self.rp
        nhat = np.sqrt(mu/(p**3))
        cots = 3.0 * nhat * self.tp
        s = np.arctan(1.0/cots)
        cottheta = (1./(np.tan(s/2.)))**0.3333
        theta = np.arctan(1./cottheta)
        tanfon2 = 2./np.tan(2.*theta)
        r = (p/2.)*(1+tanfon2**2)
        
          
        vel = np.sqrt(2.*mu/r)
        sinsqphi = p/(2.*r)
        phi = np.arcsin(np.sqrt(sinsqphi))
        f = 2.*np.arctan(tanfon2)
        xc = -r*np.cos(f)
        yc = r*np.sin(f)
        vxc = vel*np.cos(f+phi)
        vyc = -vel*np.sin(f+phi)
        xcom = self.m2 * xc/(self.m1+self.m2)
        ycom = self.m2 * yc/(self.m1+self.m2)
        vxcom = self.m2 * vxc/(self.m1 + self.m2)
        vycom = self.m2 * vyc /(self.m1+self.m2)

        self.bod1pos = np.array([[-xcom],[-ycom],[0.0]])
        self.bod1vel = np.array([[-vxcom],[-vycom],[0.0]])
        self.bod2pos = np.array([[xc-xcom],[yc-ycom],[0.0]])
        self.bod2vel = np.array([[vxc-vxcom],[vyc-vycom],[0.0]])
      



class GUI(QMainWindow):
    def __init__(self):
        super(GUI,self).__init__()
        uic.loadUi("GalCrashFinal.ui",self)
        self.show()

    
        
        #When the start button is pressed, we start seeing an animated graph in which we see how the galaxy and stars evolve over time and how they collide
        self.buttonstart.clicked.connect(self.update_plot)
        
        #The reset button clears the graph and new parameters can be selected and the simulation can be run again
        self.restartbutton.clicked.connect(self.reset)
        self.pausebutton.clicked.connect(self.pause)
        
        

        self.time = 0.00
       
        
        #graphicsView is the graph in which the galaxies are displayed
        self.graphicsView.setXRange(-20,20)
        self.graphicsView.setYRange(-20,20)
        self.graphicsView.hideAxis('bottom')
        self.graphicsView.hideAxis('left')
        
        #velplot is a plot of the relative velocities over time
        self.velplot.setXRange(0,2000)
        self.velplot.setYRange(0,1000)
        self.velplot.setLabel('bottom', 'Time [Myr]')
        self.velplot.setLabel('left', 'Relative Velocity [km/s]')
        
        #Create different colors which are used for the axes
        #mkPen outlines the item in the proposed color
        #mkBrush fills the item in the proposed color
        
        #Yellow
        self.pen1 = pg.mkPen(color = (173,171,2))
        self.brush1 = pg.mkBrush('y')
        
        #Purple/Pink
        self.pen2 = pg.mkPen(color = (155,41,163))
        
        #White
        self.pen3 = pg.mkPen(color = (255,255,255))
        
        #Used for green text
        self.pen4 = pg.mkPen(color = (85,170,0))
        #Used for red text
        self.pen5 = pg.mkPen(color = (184,6,0))
        
        #light blue
        self.pen6 = pg.mkPen(color =(28,176,217))
        self.brush6 = pg.mkBrush(color = (3,252,252))
        
        #Changing axis colors
        self.velplot.plotItem.getAxis('left').setPen(self.pen2)
        self.velplot.plotItem.getAxis('bottom').setPen(self.pen3)
        
        
        #distplot is a plot of the galaxy separation over time
        self.distplot.setXRange(0,2000)
        self.distplot.setYRange(0,200)
        #Plotted above the velocity plot so we hide the x-axis as this is common between the two
        self.distplot.hideAxis('bottom')
        self.distplot.setLabel('left','Galaxy Separation [kpc]')
        self.distplot.plotItem.getAxis('left').setPen(self.pen4)
        
        #Add labels to our graphicsView plot which give the Time, Galaxy Separation, and Relative Velocity
        self.distlabel.setStyleSheet("color: rgb(85,170,0)")
        self.distlabel.setText("Galaxy Separation:")
        self.vellabel.setStyleSheet("color: rgb(155,41,163)")
        self.vellabel.setText("Relative Velocity:")
        self.timelabel.setStyleSheet("color: white")
        self.timelabel.setText("Time:")
        
        #Creating lists for which values are appended to when update plot is active
        #These are then used to create the plots of the Galaxy Separation and Relative Velocity over Time
        # Cleared every time reset is pressed
        self.Dist = []
        self.Vel = []
        self.Time = []
        
        
        
        #Below, we configure the settings of 4 small graphs which show the impact of changing theta and phi on the orientation of the galaxy
        #Function draw calculates the orientation of an oval using theta and phi given by the user and then plots the oval
        
        self.draw()
        
        #Following commands update the image when the value of theta or phi is changed
        self.blueThetaslider.valueChanged.connect(self.draw)
        self.bluephislider.valueChanged.connect(self.draw)
        self.yellowThetaslider.valueChanged.connect(self.draw)
        self.yellowphislider.valueChanged.connect(self.draw)
        
        #Setting identical parameters for each of our graphs
        self.yellowphigraph.setXRange(-4,4)
        self.yellowphigraph.setYRange(-4,4)
        self.yellowphigraph.hideAxis('bottom')
        self.yellowphigraph.hideAxis('left')
        
        self.yellowthetagraph.setXRange(-4,4)
        self.yellowthetagraph.setYRange(-4,4)
        self.yellowthetagraph.hideAxis('bottom')
        self.yellowthetagraph.hideAxis('left')
        
        self.bluephigraph.setXRange(-4,4)
        self.bluephigraph.setYRange(-4,4)
        self.bluephigraph.hideAxis('bottom')
        self.bluephigraph.hideAxis('left')
        
        
        self.bluethetagraph.setXRange(-4,4)
        self.bluethetagraph.setYRange(-4,4)
        self.bluethetagraph.hideAxis('bottom')
        self.bluethetagraph.hideAxis('left')
        
        self.MakeGalaxy()
        self.MakeOrbit()
        
        
        self.starsnum.valueChanged.connect(self.MakeGalaxy)
        self.massratio.valueChanged.connect(self.MakeGalaxy)
        self.peribox.valueChanged.connect(self.MakeGalaxy)
        self.blueThetaslider.valueChanged.connect(self.MakeGalaxy)
        self.bluephislider.valueChanged.connect(self.MakeGalaxy)
        self.yellowThetaslider.valueChanged.connect(self.MakeGalaxy)
        self.yellowphislider.valueChanged.connect(self.MakeGalaxy)
        self.checkbighalo.stateChanged.connect(self.MakeGalaxy)
        self.bluecenterbox.stateChanged.connect(self.MakeGalaxy)
        self.yellowcenterbox.stateChanged.connect(self.MakeGalaxy)
        
        self.starsnum.valueChanged.connect(self.reset)
        self.massratio.valueChanged.connect(self.reset)
        self.peribox.valueChanged.connect(self.reset)
        self.blueThetaslider.valueChanged.connect(self.reset)
        self.bluephislider.valueChanged.connect(self.reset)
        self.yellowThetaslider.valueChanged.connect(self.reset)
        self.yellowphislider.valueChanged.connect(self.reset)
        self.checkbighalo.stateChanged.connect(self.reset)
        self.bluecenterbox.stateChanged.connect(self.reset)
        self.yellowcenterbox.stateChanged.connect(self.reset)
        
    def draw(self):
        """Function draw uses the inclination and slew angles given by the user to demonstrate how changing these angles affects the orientation of the galaxy"""
        
        #First clear graphs
        self.yellowphigraph.clear()
        self.yellowthetagraph.clear()
        self.bluephigraph.clear()
        self.bluethetagraph.clear()
        
        #Creating a regular oval with radius r in x, and radius 2r in y
        angle = np.linspace(0,2*np.pi,200)
        r = 1
        x = r * np.cos(angle)
        y = 2*r * np.sin(angle)
        z = 0
        coords = np.array([x,y,z])

        
        #In each process we take the input, convert to radians, and then use our rotation matrices to calculate new coordinates
        rtdeg = int(self.yellowThetabox.value())
        rt = rtdeg * np.pi /180
        rct = np.cos(rt)
        rst = np.sin(rt)
        Rotatert = np.array([[1,0,0],[0,rct,rst],[0,-rst,rct]])

        newrtc = Rotatert.dot(coords)
        self.yellowthetagraph.plot(newrtc[0],newrtc[1], pen = self.pen1)
        
        rpdeg = int(self.yellowphibox.value())
        rp = rpdeg * np.pi /180
        rcp = np.cos(rp)
        rsp = np.sin(rp)
        Rotaterp = np.array([[rcp,rsp,0],[-rsp,rcp,0],[0,0,1]])
        newrpc = Rotaterp.dot(coords)
        self.yellowphigraph.plot(newrpc[0],newrpc[1],pen= self.pen1)


        gtdeg = int(self.blueThetabox.value())
        gt = gtdeg * np.pi / 180
        gct = np.cos(gt)
        gst = np.sin(gt)
        Rotategt = np.array([[1,0,0],[0,gct,gst],[0,-gst,gct]])
        newgtc = Rotategt.dot(coords)
        self.bluethetagraph.plot(newgtc[0],newgtc[1],pen= self.pen6)
        
        gpdeg = int(self.bluephibox.value())
        gp = gpdeg * np.pi /180
        gcp = np.cos(gp)
        gsp = np.sin(gp)
        Rotategp = np.array([[gcp,gsp,0],[-gsp,gcp,0],[0,0,1]])
        newgpc = Rotategp.dot(coords)
        self.bluephigraph.plot(newgpc[0],newgpc[1],pen= self.pen6)
        

        
        
        
        
        
        
    def MakeGalaxy(self):
        """The function MakeGalaxy creates two galaxies based on the parameters provided in the UI"""
        #gal represents the main galaxy while comp represents the companion galaxy
        
        self.graphicsView.clear()
        self.velplot.clear()
        self.distplot.clear()
        self.graphicsView.setXRange(-20,20)
        self.graphicsView.setYRange(-20,20)
        
        #First some general parameters taken from Chris Mihos' version
        galmass = 4.8
        ahalo = 0.1
        vhalo = 1.0
        rthalo = 5.0
        galpos = np.full((3,1),0)
        galvel = np.full((3,1),0)
        diskSize = 2.5
        
        #The following values are given by the inputs the user has provided
        galtheta = int(self.blueThetabox.value())
        galphi = int(self.bluephibox.value())
        comptheta = int(self.yellowThetabox.value())
        compphi = int(self.yellowphibox.value())
        galn = int(0.5 * int(self.starsnum.value()))
        compn = int(0.5 * int(self.starsnum.value()))
        self.gal = StarGalaxy(galmass,ahalo,vhalo,rthalo,galpos,galvel,diskSize,galtheta,galphi,galn)
        self.comp = StarGalaxy(galmass,ahalo,vhalo,rthalo,galpos,galvel,diskSize,comptheta,compphi,compn)
        
        
        #If big halos is enabled, then our parameters are slightly different
        if self.checkbighalo.isChecked():
            self.gal.rthalo = 20.0
            self.gal.galmass = (self.gal.vhalo**2 * self.gal.rthalo**3)/((self.gal.ahalo+self.gal.rthalo)**2)
        else:
            self.gal.galmass = 4.8
            self.gal.rthalo = 5.0
            
        #Mass ratio is provided by the user
        massrat = float(self.massratio.value())
        
        #Using the mass ratio, we can scale the parameters of the companion galaxy according to the function scaleMass
        self.comp.scaleMass(massrat)
        
        #If big halos is enabled, then our parameters are slightly different
        if self.checkbighalo.isChecked():
            self.comp.rthalo = self.comp.rthalo*4.0
            self.comp.galmass = (self.comp.vhalo**2 * self.comp.rthalo**3)/((self.comp.ahalo+self.comp.rthalo)**2)
            
        self.MakeOrbit()
        


    def MakeOrbit(self):
        """Make orbit creates the parabolic orbit for both galaxies. It then used this orbit to determine the initial position and velocity of the galaxy centers. Finally, it initialize the stars inside both galaxies and then plots the galaxies along with their orbiting stars"""
        
        energy = 0
        eccentricity = 1
        rperi = 3.0
        
        
        #Parameter provided by user (provides the distance of closest approach)
        tperi = float(self.peribox.value())
        
        #Returns the initial position and velocities of our galaxy centers
        self.crashOrbit = Orbit(energy,rperi,tperi,eccentricity,self.gal.galmass,self.comp.galmass,self.gal.galpos,self.comp.galpos,self.gal.galvel,self.comp.galvel)
        
        #Sets the initial position and velocity of the galaxy centers
        self.gal.setPosvel(self.crashOrbit.bod1pos,self.crashOrbit.bod1vel)
        self.comp.setPosvel(self.crashOrbit.bod2pos,self.crashOrbit.bod2vel)
        
        
        #If blue centered is checked, then we alter the X and Y range of the graph so that the blue(main) galaxy
        # is always in the center of our plot
        if self.bluecenterbox.isChecked():
            self.graphicsView.setXRange(self.gal.galpos[0,0]-20,self.gal.galpos[0,0]+20)
            self.graphicsView.setYRange(self.gal.galpos[1,0]-20,self.gal.galpos[1,0]+20)
        #The same principle holds for yellow centered
        if self.yellowcenterbox.isChecked():
            self.graphicsView.setXRange(self.comp.galpos[0,0]-20,self.comp.galpos[0,0]+20)
            self.graphicsView.setYRange(self.comp.galpos[1,0]-20,self.comp.galpos[1,0]+20)
        #Plots the position of the galactic center for both galaxies
        #We give the galaxies slightly different colors to that of the stars to distinguish more easily between the two
        Gal = pg.ScatterPlotItem(brush = self.brush6, size=6)
        Gal.setData([self.gal.galpos[0,0]],[self.gal.galpos[1,0]])
        
        Comp = pg.ScatterPlotItem(brush=self.brush1, size=6)
        Comp.setData([self.comp.galpos[0,0]],[self.comp.galpos[1,0]])
        
        #Plots the points onto our plot
        self.graphicsView.addItem(Gal)
        self.graphicsView.addItem(Comp)
        
        

        
        #Initializes the parameters for the stars in both galaxies
        self.gal.InitStars()
        self.comp.InitStars()

        
        
        #Plots the surrounding stars for both galaxies
        self.graphicsView.plot((self.gal.starpos[0,:]),(self.gal.starpos[1,:]),pen = None, symbol = "t1", symbolPen = self.pen6,symbolSize = 3) 
        self.graphicsView.plot(self.comp.starpos[0,:],self.comp.starpos[1,:],pen = None, symbol = "t1", symbolPen = self.pen1, symbolSize = 3, fill=True) 
        
        #dist is the galaxy separation
        dist = 3.5 * np.linalg.norm((self.gal.galpos-self.comp.galpos))
        #We then print the galaxy separation on our graph
        self.distlabel.setText("Galaxy Separation: {:.2f} kpc".format(dist))
        
        #vel is the relative velocity of the two galaxies
        vel = 250. * np.linalg.norm((self.gal.galvel-self.comp.galvel))
        #We then print the relative velocity on our graph
        self.vellabel.setText("Relative Velocity: {:.2f} km/s".format(vel))
        
        #Print inital time of t=0
        self.timelabel.setText("Time: 0 Myr")

        

    def update_plot(self):
        """The update plot function runs the simulation. It calculates the acceleration of the galactic centers for both galaxies as well as the stars within the galaxies. These are then used to move the galaxy centers as well as the stars. The plot is then updated to show this evolution."""
        
        #We set the restartbutton and pausebutton to True here, while these are active, the plot is constantly updated
        #If the pause button is pressed, the update plot function stops running but we can still see our plots
        #If the reset button is pressed, then we clear all our graphs and reset the time back to 0
        self.restartbutton.clicked = True
        self.pausebutton.clicked = True
        
        #dtime taken from Chris Mihos' version, used to evolve the system over time
        dtime = 0.04
        if self.restartbutton.clicked:
            while self.pausebutton.clicked:
                #Clears the graphs so that the new positions can be plotted without overlapping the previous plot
                self.graphicsView.clear()
                self.distplot.clear()
                self.velplot.clear()

                #Calculating the acceleration using our acceleration function
                self.gal.galacc = self.comp.Acceleration(self.gal.galpos)
                self.comp.galacc = self.gal.Acceleration(self.comp.galpos)
                #dist is the relative distance between the two galaxies
                dist = 3.5 * np.linalg.norm((self.gal.galpos-self.comp.galpos))
                

                #If Dynamic Friction is present (given by user), then acceleration is calculated slightly differently
                if self.checkfriction.isChecked():
                    self.gal.galacc = self.gal.galacc + self.comp.DynFric(self.gal.InteriorMass(dist/3.5), self.gal.galpos,self.gal.galvel)
                    self.comp.galacc = self.comp.galacc + self.gal.DynFric(self.comp.InteriorMass(dist/3.5), self.comp.galpos,self.comp.galvel)

                comacc = ((self.gal.galmass*self.gal.galacc) + (self.comp.galmass*self.comp.galacc))/(self.gal.galmass + self.comp.galmass)
                self.gal.galacc = self.gal.galacc - comacc
                self.comp.galacc = self.comp.galacc - comacc

                self.gal.staracc = self.gal.Acceleration(self.gal.starpos) + self.comp.Acceleration(self.gal.starpos)
                self.comp.staracc = self.comp.Acceleration(self.comp.starpos) + self.gal.Acceleration(self.comp.starpos)

                #After the acceleration is calculated, we then evolve the system using the functions MoveGalaxy and MoveStars


                self.gal.MoveGalaxy(dtime)
                self.gal.MoveStars(dtime)

                self.comp.MoveGalaxy(dtime)
                self.comp.MoveStars(dtime)

                #dist and vel are the relative position and velocity between the two galaxies
                dist = 3.5 * np.linalg.norm((self.gal.galpos-self.comp.galpos))
                vel = 250. * np.linalg.norm((self.gal.galvel-self.comp.galvel))


                #Again, if blue centered is checked, then we want the graph to constantly have the blue (main) galaxy at the focus, so we constantly update the X and Y axis limits so that the galaxy appears at the center
                if self.bluecenterbox.isChecked():
                    self.graphicsView.setXRange(self.gal.galpos[0,0]-20,self.gal.galpos[0,0]+20)
                    self.graphicsView.setYRange(self.gal.galpos[1,0]-20,self.gal.galpos[1,0]+20)
                #Similarly for the yellow galaxy    
                if self.yellowcenterbox.isChecked():
                    self.graphicsView.setXRange(self.comp.galpos[0,0]-20,self.comp.galpos[0,0]+20)
                    self.graphicsView.setYRange(self.comp.galpos[1,0]-20,self.comp.galpos[1,0]+20)
                
                #Plotting the positions of the galaxies exactly the same as before
                Gal = pg.ScatterPlotItem(brush=self.brush6, symbol='o', size=6)
                Gal.setData([self.gal.galpos[0,0]],[self.gal.galpos[1,0]])
                Comp = pg.ScatterPlotItem(brush=self.brush1,symbol='o', size=6)
                Comp.setData([self.comp.galpos[0,0]],[self.comp.galpos[1,0]])
                self.graphicsView.addItem(Gal)
                self.graphicsView.addItem(Comp)
                
                #Likewise, we also plot the positions of the stars exactly the same as before
                self.graphicsView.plot((self.gal.starpos[0,:]),(self.gal.starpos[1,:]),pen = None, symbol = "t1", symbolPen = self.pen6, symbolSize = 3) 
                self.graphicsView.plot((self.comp.starpos[0,:]),(self.comp.starpos[1,:]),pen = None, symbol = "t1", symbolPen = self.pen1, symbolSize = 3)

                
                #Here we append our galaxy separation and relative velocity to the arrays defined at the start
                #These are then used to plot the galaxy separation and relative velocity as a function of time
                self.Dist.append(dist)
                self.Vel.append(vel)
                #Result is multiplied by 12 to be in Myr (taken from Chris Mihos' version)
                self.Time.append(self.time*12)
                
                #Updates the labels with the new values of each parameter
                self.distlabel.setText("Galaxy Separation: {:.2f} kpc".format(dist))
                self.vellabel.setText("Relative Velocity: {:.2f} km/s".format(vel))
                self.timelabel.setText("Time: {:.2f} Myr".format(self.time*12))
                
                #Plots the relative velocity and galaxy separation
                self.velplot.plot(self.Time,self.Vel,pen=self.pen2, width = 2, name = "Relative Velocity")
                self.distplot.plot(self.Time,self.Dist,pen=self.pen4, width = 2, name = "Galaxy Separation")

                
                #processEvents updates the plots
                QApplication.processEvents()

                #Advance the time
                self.time += dtime
        else:
            None
            
            
            
            
    def reset(self):
        """Reset clears all our visual plots and informationand resets the time to 0, so that new parameters can be given and the simulation can be run again"""
        
        self.restartbutton.clicked = False
        self.pausebutton.clicked = False
        self.Dist = []
        self.Vel = []
        self.Time = []
        self.time = 0.00
        self.distlabel.setText("Galaxy Separation:")
        self.vellabel.setText("Relative Velocity:")
        self.timelabel.setText("Time:")
        self.MakeGalaxy()

    def pause(self):
        """We create a seperate function for when pause is clicked that simply sets the pausebutton.clicked boolean to False This then means that our plot will no longer be updated until start is pressed once again"""
        self.pausebutton.clicked = False
        #When start is pressed again, we re-run our update_plot function in which we immediately set the pausebutton boolean to True again, so that our simulation will continue to run from where it left off.


        

 # here we check if there already exists an QApplication. If not we create a new one.
if not QApplication.instance(): 
    app = QApplication(sys.argv)
else:
    app = QApplication.instance()    

window = GUI() 



# here the application is started


app.exec_() #for use an in interactive shell     
