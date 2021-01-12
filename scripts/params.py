import numpy as np
def non_dimensionalize(tag,value):
	if tag == 'm_f':
		return_value = value/fixed_params['m0']
	elif tag == 'D_M1_r_h' or tag == 'D_M1_c_h' or tag == 'D_M2_h' or tag =='D_d' or tag == 'D_TNF' or tag == 'D_IL4': ## any form of a diffusion
		return_value = value*fixed_params['T']/fixed_params['L']**2
	elif tag == 'K_mM1_D' or tag == 'K_mM2_D' or tag == 'K_mM1_C':
		return_value = value*0.001/fixed_params['m0'] 
	elif tag == 'K_TNFM1_C' or tag == 'K_TNFM1_A' or tag == 'K_IL4M2_A'or tag == 'K_IL4M1_I' or tag == 'K_TGFM1_I' or tag == 'K_TNFM1_Y' or tag == 'E_TNF_h' or tag == 'E_IL4_h' :
		return_value = value*0.001/fixed_params['g0'] 
	elif tag == 'A_M1_h' or tag == 'A_M2_h' or tag == 'Y_M1_h' or tag == 'Y_M2_h' or tag == 'Y_M1_l' or tag == 'I_M1_h' or tag=='d_d' or tag == 'd_IL4' or tag == 'd_IL40' or tag == 'd_TNF' or tag == 'd_TNF0':
		return_value = value*fixed_params['T'] 
	else:
		print('non dimensionalization is not defined for {}'.format(tag))
		return_value = value
	return return_value

fixed_params = dict(
	T = 1 , ## hour
	m0 = 0.1 * 0.001, ## ng/mm3 
	L= 3.5, ## mm
	c0=10**6 * 0.001, ## cells/mm3
	g0=100 * 0.001, ## ng/mm3
	D_d = 1.23*0.0001, ## mm2/h
	D_TNF = 0.108, ## mm2/h
	D_IL4 = 0.108, ## mm2/h
	K_TNF = 0.5, ## unitless
	K_IL4 =0.5 ## unitless

)

free_params = dict(
	D_M1_r_h= [0.000036,0.0036], ## mm2/h
	D_M1_c_h= [0.000036,0.0036], ## mm2/h
	D_M2_c_h= [0.000036,0.0036], ## mm2/h
	D_M2_h= [0.000036,0.0036], ## mm2/h
	D_d= [0.000036,0.0036], ## mm2/h
	K_mM1_D = [0,0.1], ## ng/ml
	K_mM2_D = [0,0.1], ## ng/ml
	K_mM1_C = [0,0.1], ## ng/ml
	K_TNFM1_C = [0.001,1], ## ng/ml
	K_TNFM1_A = [0.001,1], ## ng/ml
	K_IL4M2_A = [0.001,1], ## ng/ml
	A_M1_h = [1/20,1/30], ## /h
	A_M2_h = [1/20,1/30], ## /h
	K_sM1_A = [0,1], ## unitless
	K_sM2_A = [0,1], ## unitless
	I_M1_h = [1/24 , 1/6],  ## /h
	K_IL4M1_I =[10,20] , ## ng/ml
	K_TGFM1_I = [5,10] , ## ng/ml
	Y_M1_h = [1/48,1/24], ## /h
	Y_M2_h = [1/48,1/24], ## /h
	Y_M1_l = [1/72,1/48], ## /h
	K_TNFM1_Y = [10,20], ## ng/ml
	E_TNF_h= [0.1,1], ## ng/ml
	E_IL4_h= [0.1,1], ## ng/ml
	d_d = [1/48,1/24], ## /h
	d_TNF0 = [1/48,1/24], ## /h
	d_IL40 = [1/48,1/24], ## /h
	d_TNF = [1/48,1/24], ## /h
	d_IL4 = [1/48,1/24], ## /h
	)



## non dimensionalize fixed params
params = {}
for tag,value in fixed_params.items():
	params.update({tag:non_dimensionalize(tag,value)})

## non dimensionalize free params
for tag,value in free_params.items():
	mean_value = np.mean(value)
	params.update({tag:non_dimensionalize(tag,mean_value)})

dependant_params = dict( # unitless
	C_TNFM1 = 2*params["K_TNFM1_C"]*params["D_M1_c_h"],
	C_mM1 = 2*params["K_mM1_C"],
	A_TNFM1 = 2*params["K_TNFM1_A"]*params["A_M1_h"],
	A_IL4M2 = 2*params["K_IL4M2_A"]*params["A_M2_h"],
	A_sM1 = 2*params["K_sM1_A"],
	A_sM2 = 2*params["K_sM2_A"],
	I_IL4M1 = 2*params["K_IL4M1_I"]*params["I_M1_h"],
	I_TGFM1 = 2*params["K_TGFM1_I"],
	)

params.update(dependant_params) 

# Conversion from Jalil's notation to mine
set_params = dict(
    # DM1
    D_mM1 = 'D_M1_r_h',
    K_mM1_D = 'K_mM1_D',
    
    # DM2
    D_mM2 = 'D_M2_h',
    K_mM2_D = 'K_mM2_D',
    
    # CM1
    C_TNFM1 = 'C_TNFM1',
    K_TNFM1_C = 'K_TNFM1_C',
    C_mM1 = 'C_mM1',
    K_mM1_C = 'K_mM1_C',
    
    # AM1
    A_TNFM1 = 'A_TNFM1',
    A_sM1 = 'A_sM1',
    K_TNFM1_A = 'K_TNFM1_A',
    K_sM1_A = 'K_sM1_A',
    
    # AM2
    A_IL4M2 = 'A_IL4M2',
    A_sM2 = 'A_sM2',
    K_IL4M2_A = 'K_IL4M2_A',
    K_sM2_A = 'K_sM2_A',
    
    # F5
    I_IL4M1 = 'I_IL4M1',
    I_TGFM1 = 'I_TGFM1',
    K_IL4M1_I = 'K_IL4M1_I',
    K_TGFM1_I = 'K_TGFM1_I',
    
    # YM1
    Y_M1_h = 'Y_M1_h',
    Y_M1_L = 'Y_M1_l',
    K_TNFM1_Y = 'K_TNFM1_Y',
    
    # YM2
    Y_M2_h = 'Y_M2_h',

    # ETNF
    E_TNF_h = 'E_TNF_h',
    K_TNF = 'K_TNF',
    
    # EIL4
    E_IL4_h = 'E_IL4_h',
    K_IL4 = 'K_IL4',
    
    # DTNF
    D_TNF_h = 'D_TNF',
    
    # DIL4
    D_IL4_h = 'D_IL4',
    
    # HTNF
    d_TNF = 'd_TNF',
    
    # HIL4
    d_IL4 = 'd_IL4',
    
    # Dd
    D_d = 'D_d',
    
    # dd
    d_d = 'd_d'
)

# Print settings to paste them derictly to settings file
print('\nPaste the following directly to settings file, subsection "Physical constants":\n')
for key in set_params:
    print("set {} = {}".format(key, params[set_params[key]]))
    