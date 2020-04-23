##
# Newton-Euler Recursive Method Calculation
# for Robotics Dynamic Implentation
#
# @author Rafael Lima
##

from sympy import *

# Generic DH Parameters
thetai = Symbol('theta_i');
ALPHAi = Symbol('ALPHA_i');
di = Symbol('d_i');
ai = Symbol('a_i');

# Generic DH Parameters Transformation Matrix
A_i = Matrix([[cos(thetai),-cos(ALPHAi)*sin(thetai),sin(ALPHAi)*sin(thetai),ai*cos(thetai)],
               [sin(thetai),cos(ALPHAi)*cos(thetai),-sin(ALPHAi)*cos(thetai),ai*sin(thetai)],
               [0,sin(ALPHAi),cos(ALPHAi),di],
               [0,0,0,1]]);

# Generic Rotation Matrix
R_i = Matrix([[cos(thetai),-cos(ALPHAi)*sin(thetai),sin(ALPHAi)*sin(thetai)],
                [sin(thetai),cos(ALPHAi)*sin(thetai),-sin(ALPHAi)*cos(thetai)],
                [0,sin(ALPHAi),cos(ALPHAi)]]);

# Generic Position Matrix
P_i = Matrix([ai*cos(thetai),ai*sin(thetai),di]);

##
# 2 Links Robots - Kinematics
##

# Link 1 - DH Parameters:
theta1 = Symbol('theta_1');
ALPHA1 = 0;
d1 = 0;
a1 = Symbol('a_1');

# Link 1 - Transformation Matrix:
A_01 = A_i.subs([(thetai,theta1),
                 (ALPHAi,ALPHA1),
                 (di,d1),
                 (ai,a1)]);

# Link 1 - Rotation Matrix
R_01 = R_i.subs([(thetai,theta1),
                 (ALPHAi,ALPHA1),
                 (di,d1),
                 (ai,a1)]);

R_10 = R_01.T # Transpose

# Link 1 - Position Matrix
P_01 = P_i.subs([(thetai,theta1),
                 (ALPHAi,ALPHA1),
                 (di,d1),
                 (ai,a1)]);

# Link 1 - Position Vector Center of Gravity 
Xc1 = Symbol('X_c1');
Yc1 = Symbol('Y_c1');
Zc1 = Symbol('Z_c1');

P_c1 = Matrix([Xc1,Yc1,Zc1]);

# Link 2 - DH Parameters:
theta2 = Symbol('theta_2');
ALPHA2 = 0;
d2 = 0;
a2 = Symbol('a_2');

# Link 2 - Transformation Matrix:
A_12 = A_i.subs([(thetai,theta2),
                 (ALPHAi,ALPHA2),
                 (di,d2),
                 (ai,a2)]);

# Link 2 - Rotation Matrix
R_12 = R_i.subs([(thetai,theta2),
                 (ALPHAi,ALPHA2),
                 (di,d2),
                 (ai,a2)]);

R_21 = R_12.T; # Transpose

# Link 2 - Position Vector
P_12 = P_i.subs([(thetai,theta2),
                 (ALPHAi,ALPHA2),
                 (di,d2),
                 (ai,a2)]);

# Link 2 - Position Vector Center of Gravity 
Xc2 = Symbol('X_c2');
Yc2 = Symbol('Y_c2');
Zc2 = Symbol('Z_c2');

P_c2 = Matrix([Xc2,Yc2,Zc2]);

# Link 1 - Angular Velocity
omega1 = Symbol('omega_1');
omega_01 = Matrix([0,0,omega1]);

# Link2 - Angular Velocity
omega2 = Symbol('omega_2');
omega_12 = Matrix([0,0,omega2]);

##
# Foward Computation for Velocity and Acceleration
# for Revolute Joints
##

# Angular and Linear Velocity for revolute Joints

omega_101 = R_10*omega_01;

v_101 = omega_101.cross(R_10*P_01);

v_101.simplify();

omega_202 = R_21*(omega_101 + omega_12);

v_202 = R_21*v_101 + omega_202.cross(R_21*P_12);


# Angular Acceleration Propagation

alpha1 = Symbol('alpha_1');
alpha_01 = Matrix([0,0,alpha1]);

alpha2 = Symbol('alpha_2');
alpha_12 = Matrix([0,0,alpha2]);

alpha_101 = R_10*alpha_01;
alpha_202 = R_21*(alpha_101 + omega_101.cross(omega_12) + alpha_12);

# Linear Acceleration Propagation:
a_101 = R_10*(alpha_01.cross(P_01) + omega_01.cross(omega_01.cross(P_01)));

a_202 = R_21*(a_101+ alpha_12.cross(P_12) + omega_12.cross(omega_12.cross(P_12)) + 2*omega_101.cross(omega_12.cross(P_12)) + alpha_101.cross(P_12) + omega_101.cross(omega_101.cross(P_12)));

# Unable to simplify a_101 , a_202

# Radial Distance to the center of mass of each link
a1 = Symbol('a_1');
a2 = Symbol('a_2');

P_101 = Matrix([a1,d1*sin(theta1),d1*cos(theta1)]);
P_212 = Matrix([a2,d2*sin(theta2),d2*cos(theta2)]);

P_c1 = Matrix([-a1/2,0,0]);
P_c2 = Matrix([-a2/2,0,0]);

# Acceleration at the center of mass of each link:

a_c101 = a_101 + alpha_101.cross(P_c1) + omega_101.cross(omega_101.cross(P_c1));

a_c202 = a_202 + alpha_202.cross(P_c2) + omega_202.cross(omega_202.cross(P_c2));

##
# Backward Computation of forces and moments
##

# Determination of the gravity vectors

# g= 9.81
g = Symbol('g')

GR_0 = Matrix([0,-g,0]);
GR_1 = R_10*GR_0;

R_20 = (R_01*R_12).T; # TODO Calculate R_20;

GR_2 = R_20*GR_0;


## Determine Inercia Matrix

m1 = Symbol('m_1');
m2 = Symbol('m_2');

II1 = 1.0/12*(m1*a1**2)*Matrix([[0,0,0],[0,1,0],[0,0,1]]);
II2 = 1.0/12*(m2*a2**2)*Matrix([[0,0,0],[0,1,0],[0,0,1]]);

# Inercia Force and Moments
f_22tool = Matrix([0,0,0]);
n_22tool = Matrix([0,0,0]);

Z = Matrix([0,0,1]);

f_22 = m2*a_c202;
n_22 = II2*alpha_202 - omega_202.cross(II2*omega_202);

f_212 = f_22tool - f_22 - m2*GR_2;
f_112 = R_12*f_212;

n_212 = n_22tool + (P_212+P_c2).cross(f_212) - P_c2.cross(f_22tool) - n_22;
n_112 = R_12*n_212;

n_11 = II1*alpha_101 - omega_101.cross(II1*omega_101);
f_11 = m1*a_c101;

f_101 = f_112 - f_11 - m1*GR_1;
n_101 = n_112+(P_101+P_c1).cross(f_101) - P_c1.cross(f_112) - n_11;

n_001 = R_01*n_101;

##
# Dynamic Equation
##
