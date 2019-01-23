import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from scipy.fftpack import fft, ifft, diff
from scipy.integrate import odeint, complex_ode

alpha = 1
beta = 1
grid = 1
dom = 200
t0 = 0
t1 = 400
dt = 1

def dAdt(t, A):

	L = len(A)
	Re_A_xx = diff(A.real, order=2, period=L)	
	Im_A_xx = diff(A.imag, order=2, period=L)	
	A_xx = Re_A_xx + 1j*Im_A_xx
	abs_A2 = abs(A)**2

	A_t = (1 + 1j*alpha)*A_xx + A - (1 + 1j*beta)*abs_A2*A
	A_t_ = 1j*(A_xx + 2*abs(A)*A)

	return A_t


#sol = np.array([], dtype = np.complex64)
sol = []
t   = np.array([], dtype = np.complex64)
x = np.arange(-dom,dom,grid)
L = len(x)

mu = 0
sigma = 0.1

#A_0 = 0.2*np.sin(2*np.pi*x) + 0.1*1j*np.cos(2*np.pi*x+np.pi/6) + 0.1*np.random.random()
#A_0 = np.random.normal(mu, sigma, dom*2)
#A_0 = 1/np.cos((2*np.pi*x/L)**2) + 0.8/np.cos((2*np.pi*x/L)**2) + 0.01*np.random.random()
#A_0 = np.cos(2*np.pi*x/L) + 0.1*np.random.normal(mu, sigma, dom*2)
A_0 = np.sqrt(1 - (20*np.pi/L)**2)*np.exp(1j*20*np.pi*x/L) + 0.1*np.random.normal(mu, sigma, dom*2)
L = len(A_0)

r = complex_ode(dAdt)
r.set_initial_value(A_0, t0)

while r.successful() and r.t < t1:
	t = np.append(t,r.t+dt)
	sol.append(r.integrate(r.t+dt))

sol = np.array(sol)

#abs_sol = sol[:,0:L/2]**2 + sol[:,L/2:L]**2
#abs_sol = np.array([abs_sol,abs_sol])
#abs_sol = abs_A2.flatten()

X,T = np.meshgrid(x,t)

fig = plt.figure(figsize=(15,6))

# `ax` is a 3D-aware axis instance, because of the projection='3d' keyword argument to add_subplot
ax = fig.add_subplot(1, 2, 1, projection='3d')
# surface_plot with color grading and color bar
p = ax.plot_surface(X, T, abs(sol), rstride=1, cstride=1, cmap=plt.get_cmap('jet'), linewidth=0, antialiased=False)
#cb = fig.colorbar(p, shrink=0.5)

#plt.subplot(121)
ax = fig.add_subplot(1, 2, 2)
ax.imshow(abs(sol[::-1,:]), aspect='equal', cmap=plt.get_cmap('jet'))
#plt.imshow(sol[::-1,:].real)
fig.colorbar(p, shrink=0.5)

#plt.subplot(122)
#plt.plot(x,abs(A_0))
#plt.plot(x,abs(sol[-1,:]))
#plt.plot(x,abs(sol[2,:]))


plt.show()