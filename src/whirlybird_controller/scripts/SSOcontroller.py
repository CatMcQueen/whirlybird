#!/usr/bin/env python
# license removed for brevity


# This file is a basic structure to write a controller that
# communicates with ROS. It will be the students responsibility
# tune the gains and fill in the missing information

# As an example this file contains PID gains, but other
# controllers use different types of gains so the class
# will need to be modified to accomodate those changes

import rospy
import time
import numpy as np
import control as cnt
from whirlybird_msgs.msg import Command
from whirlybird_msgs.msg import Whirlybird
from std_msgs.msg import Float32


class Controller():

    def __init__(self):

        # get parameters
        try:
            param_namespace = '/whirlybird'
            self.param = rospy.get_param(param_namespace)
        except KeyError:
            rospy.logfatal('Parameters not set in ~/whirlybird namespace')
            rospy.signal_shutdown('Parameters not set')
	    print("Key Error Exception thrown!")

        g = self.param['g']
        l1 = self.param['l1']
        l2 = self.param['l2']
        m1 = self.param['m1']
        m2 = self.param['m2']
        d = self.param['d']
        h = self.param['h']
        r = self.param['r']
        Jx = self.param['Jx']
        Jy = self.param['Jy']
        Jz = self.param['Jz']
        km = self.param['km']
	
	#sigma for dirty derivative
	self.sigma = 0.05

        # Roll Gains
        self.P_phi_ = 0.0
        self.I_phi_ = 0.0
        self.D_phi_ = 0.0
        self.Int_phi = 0.0
        self.prev_phi = 0.0

	self.phi_dot = 0.0

        # Pitch Gains
        self.theta_r = 0.0
        self.P_theta_ = 0.0
        self.I_theta_ = 0.0
        self.D_theta_ = 0.0
        self.prev_theta = 0.0
        self.Int_theta = 0.0

	self.theta_dot = 0.0
	self.thetae = 0.0
	self.theta_observer = 0.0

	self.prev_theta_observer = 0.0

        # Yaw Gains
        self.psi_r = 0.0
        self.P_psi_ = 0.0
        self.I_psi_ = 0.0
        self.D_psi_ = 0.0
        self.prev_psi = 0.0
        self.Int_psi = 0.0
	self.psi_dot = 0.0
	self.psi_observer = 0.0

	self.prev_psi_observer = 0.0

	self.prev_time = rospy.Time.now()

	theta_init = 0.0
	self.dirtyd_theta = 0.0
	self.prev_dirtyd_theta = 0.0
	self.dirtyd_psi = 0.0
	self.prev_dirtyd_psi = 0.0
        self.Fe = (m1*l1 - m2*l2)*g*np.cos(theta_init)/l1

	self.Klat = 0.0
	self.Klon = 0.0
	self.K1lon = 0.0
	self.K1lat = 0.0
	self.lat_integrator = 0.0

	self.kilon = 0.0
	self.kilat = 0.0
	self.lon_integrator = 0.0

	self.x_lon_hat = np.matrix([[0.0],
				    [0.0]])

	self.x_lat_hat = np.matrix([[0.0],
				    [0.0],
				    [0.0],
				    [0.0]])


	zeta = .7
	tr = 1.0 #1.4
	numer = 2.2
	wn = numer/tr
	b0 = 1.152

	tr_phi = .3
	tr_psi = 5*tr_phi
	zeta_lat = .8
	numer_lat = 2.2
	wn_phi = numer_lat/tr_phi
	wn_psi = numer_lat/tr_psi
	bpsi = l2 * self.Fe / (m1*l1**2 + m2*l2**2 + Jz)

	#observer
	zeta_o_lon = .707
	tr_o_lon = tr
	wn_o_lon = 2.2 / (tr_o_lon / 5)

	zeta_o_lat = 0.707
	tr_phi_o_lat = tr_phi
	wn_phi_o_lat = 2.2 / (tr_phi_o_lat / 5)

	tr_psi_o_lat = tr_psi
	wn_psi_o_lat = 2.2 / (tr_psi_o_lat / 5)

        self.int_pol_lon = -wn/2.0
        self.int_pol_lat = -wn_psi/2.0

	#longitude Matrices
	Alon = np.matrix([[0.0, 1.0],
			  [(m1*l1-m2*l2)*g*np.sin(self.thetae)/(m1*(l1**2)+m2*(l2**2)+Jy), 0.0]])

	Blon = np.matrix([[0.0],
			  [l1/(m1*(l1**2)+m2*(l2**2)+Jy)]])

	Clon = np.matrix([[1.0, 0.0]])

	A1lon = np.matrix([[0.0, 1.0, 0.0],
			  [(m1*l1-m2*l2)*g*np.sin(self.thetae)/(m1*(l1**2)+m2*(l2**2)+Jy), 0.0, 0.0],
			  [-1.0, 0.0, 0.0]])

	B1lon = Blon = np.matrix([[0.0],
			  [l1/(m1*(l1**2)+m2*(l2**2)+Jy)],
			  [0.0]])
	
	des_char_poly_lon = np.convolve([1, 2*zeta*wn, wn**2], np.poly(self.int_pol_lon))
	des_poles_lon = np.roots(des_char_poly_lon)

	if np.linalg.matrix_rank(cnt.ctrb(A1lon, B1lon)) != 3:
		print("The longitudinal system is not controllable")
	else:
		self.K1lon = cnt.place(A1lon, B1lon, des_poles_lon)
		self.Klon = np.matrix([self.K1lon.item(0), self.K1lon.item(1)])
		self.kilon = self.K1lon.item(2)
		print("Klon: ", self.Klon, "kilon: ", self.kilon)

	lon_des_obs_char_poly = [1, 2*zeta_o_lon*wn_o_lon, wn_o_lon**2]
	lon_des_obs_poles= np.roots(lon_des_obs_char_poly)

	if np.linalg.matrix_rank(cnt.ctrb(Alon.T, Clon.T)) != 2:
		print("The systems is not observable")
	else:
		self.Llon = signal.place_poles(Alon.T, Clon.T, lon_des_obs_poles).gain_matrix.T

	#latitude Matrices

	Alat = np.matrix([[0.0, 0.0, 1.0, 0.0],
			  [0.0, 0.0, 0.0, 1.0],
			  [0.0, 0.0, 0.0, 0.0],
			  [l1*self.Fe/(m1*(l1**2)+m2*(l2**2)+Jz), 0.0, 0.0, 0.0]])

	Blat = np.matrix([[0.0],
			  [0.0],
			  [1/Jx],
			  [0.0]])

	Clat = np.matrix([[1.0, 0.0, 0.0, 0.0],
			  [0.0, 1.0, 0.0, 0.0]])

	A1lat = np.matrix([[0.0, 0.0, 1.0, 0.0, 0.0],
			  [0.0, 0.0, 0.0, 1.0, 0.0],
			  [0.0, 0.0, 0.0, 0.0, 0.0],
			  [l1*self.Fe/(m1*(l1**2)+m2*(l2**2)+Jz), 0.0, 0.0, 0.0, 0.0],
			  [0.0, -1.0, 0.0, 0.0, 0.0]])

	B1lat = Blat = np.matrix([[0.0],
			  [0.0],
			  [1/Jx],
			  [0.0],
			  [0.0]])

	des_char_poly_lat = np.convolve(np.convolve([1, 2*zeta*wn_psi, wn_psi**2],
					[1, 2*zeta*wn_phi, wn_phi**2]), np.poly(self.int_pol_lat))
	
	des_poles_lat = np.roots(des_char_poly_lat)

	if np.linalg.matrix_rank(cnt.ctrb(A1lat, B1lat)) != 5:
		print("The lateral system is not controllable")
	else:
		self.K1lat = cnt.place(A1lat, B1lat, des_poles_lat)
		self.Klat = np.matrix([self.K1lat.item(0), self.K1lat.item(1), self.K1lat.item(2), self.K1lat.item(3)])
		self.kilat = self.K1lat.item(4)

	lat_des_obs_char_poly = np.convolve([1, 2*zeta_o_lat*wn_psi_o_lat, wn_psi_o_lat**2],
					    [1, 2*zeta_o_lat*wn_phi_o_lat, wn_phi_o_lat**2])
	lat_des_obs_poles= np.roots(lat_des_obs_char_poly)

	if np.linalg.matrix_rank(cnt.ctrb(Alat.T, Clat.T)) != 2:
		print("The systems is not observable")
	else:
		self.Llat = signal.place_poles(Alat.T, Clat.T, lat_des_obs_poles).gain_matrix.T
		
        self.command_sub_ = rospy.Subscriber('whirlybird', Whirlybird, self.whirlybirdCallback, queue_size=5)
        self.psi_r_sub_ = rospy.Subscriber('psi_r', Float32, self.psiRCallback, queue_size=5)
        self.theta_r_sub_ = rospy.Subscriber('theta_r', Float32, self.thetaRCallback, queue_size=5)
        self.command_pub_ = rospy.Publisher('command', Command, queue_size=5)
	#observer stuff
	self.psi_o_sub_ =  rospy.Subscriber('psi_o', Float32, self.psiOCallback, queue_size = 5)
	self.theta_o_sub_ = rospy.Subscriber('theta_o', Float32, self.thetaOCallback, queue_size = 5)

	self.c = 0
        while not rospy.is_shutdown():
            # wait for new messages and call the callback when they arrive
            rospy.spin()
    
    def psiOCallback(self, msg):
	self.psi_observer = msg.data

    def thetaOCallback(self, msg):
	self.theta_observer = msg.data

'''
    def updateLonObserver(self, y_m):
	N = 5
	y = np.matrix([[y_m[?]]])

    def updateLatObserver():
'''

    def thetaRCallback(self, msg):
        self.theta_r = msg.data

    def psiRCallback(self, msg):
        self.psi_r = msg.data

    def whirlybirdCallback(self, msg):
        g = self.param['g']
        l1 = self.param['l1']
        l2 = self.param['l2']
        m1 = self.param['m1']
        m2 = self.param['m2']
        d = self.param['d']
        h = self.param['h']
        r = self.param['r']
        Jx = self.param['Jx']
        Jy = self.param['Jy']
        Jz = self.param['Jz']
        km = self.param['km']

        phi = msg.roll
        theta = msg.pitch
        psi = msg.yaw
	theta_observer = 
	psi_observer =

        # Calculate dt (This is variable)
        now = rospy.Time.now()
        dt = (now-self.prev_time).to_sec()
	beta = -(2.0*self.sigma - dt)/(2.0*self.sigma + dt)
        self.prev_time = now
        ##################################
        # Controller Implemented here

	#longitudinal

		

	lon_error = -theta_observer + self.theta_r
	lon_prev_error = -self.prev_theta_observer + self.theta_r
	lon_error_dot = lon_error - lon_prev_error
	
	self.lon_integrator = self.lon_integrator + ((dt/2.0) * (lon_error + lon_prev_error))
	
	self.Fe = (m1*l1 - m2*l2)*g*np.cos(theta)/l1
	
	self.theta_dot = beta*self.theta_dot + (1-beta)*((theta - self.prev_theta)/dt)

	xlon = np.matrix([[theta], [self.theta_dot]])

	Ftilde = -self.Klon*xlon - self.kilon*self.lon_integrator
	F = self.Fe + Ftilde

	#lateral
	lat_error = -psi_observer + self.psi_r
	lat_prev_error = -self.prev_psi + self.psi_r
	lat_error_dot = (lat_error - lat_prev_error)/dt
	

	self.lat_integrator = self.lat_integrator + dt/2.0 * (lat_error + lat_prev_error)

	self.phi_dot = beta*self.phi_dot + (1-beta)*((phi - self.prev_phi)/dt)
	self.psi_dot = beta*self.psi_dot + (1-beta)*((psi - self.prev_psi)/dt)
	
	xlat = np.matrix([[phi], [psi], [self.phi_dot], [self.psi_dot]])

	tau = -self.Klat*xlat - self.kilat*self.lat_integrator #+ 0.2 # tau_e = 0

	# this should be in terms of tau and F.
	T1 = tau[0] / d

	left_force, right_force = F/2+T1, F/2-T1

	#print ("left_force:", left_force)
	self.prev_theta = theta
	self.prev_psi = psi
	self.prev_phi = phi
	self.prev_theta_observer = theta_observer

        ##################################

        # Scale Output
        l_out = left_force/km
        if(l_out < 0):
            l_out = 0
        elif(l_out > 0.7):
            l_out = 0.7 

        r_out = right_force/km
        if(r_out < 0):
            r_out = 0
        elif(r_out > 0.7):
            r_out = 0.7

	#print("l_out: ", l_out, "r_out: ", r_out)

        # Pack up and send command
        command = Command()
        command.left_motor = l_out
        command.right_motor = r_out
        self.command_pub_.publish(command)


if __name__ == '__main__':
    rospy.init_node('controller', anonymous=True)
    controller = Controller()
    try:
        controller = Controller()
    except:
        rospy.ROSInterruptException
    pass
