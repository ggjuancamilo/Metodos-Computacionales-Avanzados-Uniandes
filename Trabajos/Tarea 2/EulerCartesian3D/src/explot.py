import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os.path

if os.path.isfile('data1'):
	data1 = np.loadtxt('data1')
	z1 = data1[:,2]-64
	y1 = data1[:,1]-64
	x1 = data1[:,0]-64
	rho1 = data1[:,3]	
	for i in range(np.size(x1)):
		if (np.isnan(x1[i])or(x1[i]==inf)):
			x1[i] = 1.0
		if (np.isnan(y1[i])or(y1[i]==inf)):
			y1[i] = 1.0
		if (np.isnan(z1[i])or(z1[i]==inf)):
			z1[i] = 1.0
		if (np.isnan(rho1[i])or(rho1[i]==inf)):
			rho1[i] = 1.0
			
	r = np.linspace(0.0, 126, 64)

	rhon1 = np.array([])
	for i in range(np.size(1, r)):
		temp1 = np.empty()
		for j in range(np.size(x1)):
			if ((((x1[j]**2.0)+(y1[j]**2.0)+(z1[j]**2.0))**0.5)<=r[i]) and ((((x1[j]**2.0)+(y1[j]**2.0)+(z1[j]**2.0))**0.5)>r[i-1]):
				temp1 = np.append(temp1, rho1[j])
		rhon1 = np.append(rhon1, temp1)
	
	for i in range(rhon1):
		rhon1[i] = np.mean(rhon1[i])
	rhon1[0] = rhon1[1]
	plt.plot(r, rhon1)
	plt.savefig('sedov1.pdf')

else:
	xx = np.linspace(0,126,100)
	plt.plot(xx, np.sin(xx) )
	plt.savefig('sedov1.pdf')		

if os.path.isfile('data2'):
	data2 = np.loadtxt('data2')
	x2 = data2[:,0]-64
	y2 = data2[:,1]-64
	z2 = data2[:,2]-64
	rho2 = data2[:,3]
	for i in range(np.size(x1)):
		if (np.isnan(x2[i])or(x2[i]==inf)):
			x2[i] = 1.0
		if (np.isnan(y2[i])or(y2[i]==inf)):
			y2[i] = 1.0
		if (np.isnan(z2[i])or(z2[i]==inf)):
			z2[i] = 1.0
		if (np.isnan(rho2[i])or(rho2[i]==inf)):
			rho2[i] = 1.0
	rhon2 = np.array([])

	for i in range(np.size(1, r)):
		temp1 = np.empty()
		temp2 = np.empty()
		temp3 = np.empty()
		for j in range(np.size(x1)):
			if ((((x2[j]**2.0)+(y2[j]**2.0)+(z2[j]**2.0))**0.5)<=r[i]) and ((((x2[j]**2.0)+(y2[j]**2.0)+(z2[j]**2.0))**0.5)>r[i]):
				temp2 = np.append(temp2, rho2[j])
			rhon2 = np.append(rhon2, temp2)
	for i in range(rhon1):
		rhon2[i] = np.mean(rhon2[i])
	rhon2[0] = rhon2[1]
	plt.plot(r, rhon2)
	plt.savefig('sedov2.pdf')
	
else:
	xx = np.linspace(0,126,100)
	plt.plot(xx, np.sin(xx) )
	plt.savefig('sedov2.pdf')
	
if os.path.isfile('data3'):
	data3 = np.loadtxt('data3')
	x3 = data3[:,0]-64
	y3 = data3[:,1]-64
	z3 = data3[:,2]-64
	rho3 = data3[:,3]
	for i in range(np.size(x1)):
		if (np.isnan(x3[i])or(x3[i]==inf)):
			x3[i] = 1.0
		if (np.isnan(y3[i])or(y3[i]==inf)):
			y3[i] = 1.0
		if (np.isnan(y3[i])or(y3[i]==inf)):
			y3[i] = 1.0
		if (np.isnan(z3[i])or(z3[i]==inf)):
			z3[i] = 1.0	
		if (np.isnan(rho3[i])or(rho3[i]==inf)):
			rho3[i] = 1.0
		rhon3 = np.array([])
	
	for i in range(np.size(1, r)):
		temp1 = np.empty()
		temp2 = np.empty()
		temp3 = np.empty()
		for j in range(np.size(x1)):		
			if  ((((x3[j]**2.0)+(y3[j]**2.0)+(z3[j]**2.0))**0.5)<=r[i]) and ((((x3[j]**2.0)+(y3[j]**2.0)+(z3[j]**2.0))**0.5)>r[i]):
				temp3 = np.append(temp3, rho3[j])
			rhon3 = np.append(rhon3, temp3)
	for i in range(rhon1):
		rhon1[i] = np.mean(rhon1[i])
		rhon2[i] = np.mean(rhon2[i])
		rhon3[i] = np.mean(rhon3[i])
	rhon3[0] = rhon3[1]
	plt.plot(r, rhon3)
	plt.savefig('sedov3.pdf')
	
else:
	xx = np.linspace(0,126,100)
	plt.plot(xx, np.sin(xx) )
	plt.savefig('sedov3.pdf')
