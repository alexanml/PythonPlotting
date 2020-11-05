"""
Module for the verification of entropy production calculated in cott vs analytically.
Author: AML
Oct/Nov 2020
"""
import os.path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns  # makes plots look extra pretty

import knuplot
#import pycott # Uncomment this if you have a symbolic link to pycott
import pycott_not_link as pycott

sns.set_context("paper", font_scale=1.4) 



def get_cott_sigma_vs_t_at_x(tec_file, min_timestep, max_timestep, noneq_relax_const_K=0, noneq_relax_const_H=0, friction=False):
    
    x_from_tec = tec_file.ts[0]['x']
    print("x", x_from_tec)
    dx = x_from_tec[4]-x_from_tec[3]
    
    #icell_a = (np.abs(x_from_tec - x_a)).argmin()
    #icell_b = (np.abs(x_from_tec - x_b)).argmin()
    
    #actual_x_a = x_from_tec[icell_a]
    #actual_x_b = x_from_tec[icell_b]
    #print("Wanted position, x_a:", x_a)
    #print("Closest position found in tec file, x_tec:", actual_x_a)
    #print("Wanted position, x_b:", x_b)
    #print("Closest position found in tec file, x_tec:", actual_x_b)
    
    timesteps = tec_file.nzones #number of time steps in tec-file
    tot_dSirrdt_num = np.zeros(max_timestep-min_timestep)
    tot_dSirrdt_calc = np.zeros(max_timestep-min_timestep)

    # plot IC
    p = tec.ts[min_timestep-1]['p_mix']
    plt.plot(x_from_tec,p)
    plt.xlabel("x")
    plt.ylabel(r"p")
    plt.show()

    u1 = tec.ts[min_timestep-1]['u_1']
    plt.plot(x_from_tec,u1)
    plt.xlabel("x")
    plt.ylabel(r"u")
    plt.show()

    T = tec.ts[min_timestep-1]['Temp_1']
    plt.plot(x_from_tec,T)
    plt.xlabel("x")
    plt.ylabel(r"T")
    plt.show()

    rho1 = tec.ts[min_timestep-1]["`r_1"] 
    alpha1 = tec.ts[min_timestep-1]["`a_1"]
    rho2 = tec.ts[min_timestep-1]["`r_2"] 
    alpha2 =  tec.ts[min_timestep-1]["`a_2"]
    rhotot = rho1*alpha1+rho2*alpha2
    plt.plot(x_from_tec,rhotot)
    plt.xlabel("x")
    plt.ylabel(r"$\rho$")
    plt.show()

    # plot end condition
    p = tec.ts[max_timestep-1]['p_mix']
    plt.plot(x_from_tec,p)
    plt.xlabel("x")
    plt.ylabel(r"p")
    plt.show()

    u1 = tec.ts[max_timestep-1]['u_1']
    plt.plot(x_from_tec,u1)
    plt.xlabel("x")
    plt.ylabel(r"u")
    plt.show()

    T = tec.ts[max_timestep-1]['Temp_1']
    plt.plot(x_from_tec,T)
    plt.xlabel("x")
    plt.ylabel(r"T")
    plt.show()

    rho1 = tec.ts[max_timestep-1]["`r_1"] 
    alpha1 = tec.ts[max_timestep-1]["`a_1"]
    rho2 = tec.ts[max_timestep-1]["`r_2"] 
    alpha2 =  tec.ts[max_timestep-1]["`a_2"]
    rhotot = rho1*alpha1+rho2*alpha2
    plt.plot(x_from_tec,rhotot)
    plt.xlabel("x")
    plt.ylabel(r"$\rho$")
    plt.show()

    # this is just to get the times in the simulation:
    t, irrelevantvar     =  pycott.cott_utils.get_var_vs_t("Temp_2", 1, tec_file)
    if (timesteps > 2):
            for i in range(min_timestep,max_timestep-1):
                # values at t_{i-1}
                x, s_dwn = pycott.cott_utils.get_var_vs_x("entrmix", t[i-1], tec_file)  # = s [J/kg]
                x, u_dwn = pycott.cott_utils.get_var_vs_x("u_1", t[i-1], tec_file)      # u1 = u2 = u
                x, rho1_dwn = pycott.cott_utils.get_var_vs_x("`r_1", t[i-1], tec_file) 
                x, alpha1_dwn = pycott.cott_utils.get_var_vs_x("`a_1", t[i-1], tec_file)
                x, rho2_dwn = pycott.cott_utils.get_var_vs_x("`r_2", t[i-1], tec_file)
                x, alpha2_dwn = pycott.cott_utils.get_var_vs_x("`a_2", t[i-1], tec_file)
                x, Saru_dwn       = pycott.cott_utils.get_var_vs_x("`S`a`ru", t[i-1], tec_file)
                x, p_dwn      = pycott.cott_utils.get_var_vs_x("p_mix", t[i-1], tec_file)
                x, gasfrac_dwn = pycott.cott_utils.get_var_vs_x("gasfrac", t[i-1], tec_file)
                x, mu1_dwn   = pycott.cott_utils.get_var_vs_x("`m_1", t[i-1], tec_file)
                x, mu2_dwn   = pycott.cott_utils.get_var_vs_x("`m_2", t[i-1], tec_file)
                x, mu_int_dwn =  pycott.cott_utils.get_var_vs_x("`m_mix", t[i-1], tec_file)
                x, T1_dwn     =  pycott.cott_utils.get_var_vs_x("Temp_1", t[i-1], tec_file)
                x, T2_dwn     =  pycott.cott_utils.get_var_vs_x("Temp_2", t[i-1], tec_file)
                # this is hopefully heat from the wall # no heat in this case (noneq vp)
                #x, Q_w = pycott.cott_utils.get_var_vs_x("heat", t[i], tec_file)
                # this is hopefully temperature in the wall
                #x, T_w  = pycott.cott_utils.get_var_vs_x("T_w", t[i], tec_file)  

                if friction:
                    t, fric_dwn   = pycott.cott_utils.get_var_vs_t("fric", t[i-1], tec_file)
            
                rho_dwn =  rho1_dwn*alpha1_dwn + rho2_dwn*alpha2_dwn
                rhos_tot_dwn = np.sum(s_dwn*rho_dwn)
                
                # values at t_i
                x, s = pycott.cott_utils.get_var_vs_x("entrmix", t[i], tec_file)  # = s [J/kg]
                x, u = pycott.cott_utils.get_var_vs_x("u_1", t[i], tec_file)      # u1 = u2 = u
                x, rho1 = pycott.cott_utils.get_var_vs_x("`r_1", t[i], tec_file) 
                x, alpha1 = pycott.cott_utils.get_var_vs_x("`a_1", t[i], tec_file)
                x, rho2 = pycott.cott_utils.get_var_vs_x("`r_2", t[i], tec_file)
                x, alpha2 = pycott.cott_utils.get_var_vs_x("`a_2", t[i], tec_file)
                x, Saru   = pycott.cott_utils.get_var_vs_x("`S`a`ru", t[i], tec_file)
                x, p      = pycott.cott_utils.get_var_vs_x("p_mix", t[i], tec_file)
                x, gasfrac = pycott.cott_utils.get_var_vs_x("gasfrac", t[i], tec_file)
                x, mu1   = pycott.cott_utils.get_var_vs_x("`m_1", t[i], tec_file)
                x, mu2   = pycott.cott_utils.get_var_vs_x("`m_2", t[i], tec_file)
                x, mu_int =  pycott.cott_utils.get_var_vs_x("`m_mix", t[i], tec_file)
                x, T1     =  pycott.cott_utils.get_var_vs_x("Temp_1", t[i], tec_file)
                x, T2     =  pycott.cott_utils.get_var_vs_x("Temp_2", t[i], tec_file)
                # this is hopefully heat from the wall # no heat in this case (noneq vp)
                #x, Q_w = pycott.cott_utils.get_var_vs_x("heat", t[i], tec_file)
                # this is hopefully temperature in the wall
                #x, T_w  = pycott.cott_utils.get_var_vs_x("T_w", t[i], tec_file)  

                if friction:
                    t, fric   = pycott.cott_utils.get_var_vs_t("fric", t[i], tec_file)
            
                rho = rho1*alpha1 + rho2*alpha2
                rhos_tot = np.sum(s*rho)

                # values at t_{i+1}
                x, s_up = pycott.cott_utils.get_var_vs_x("entrmix", t[i+1], tec_file)  # = s [J/kg]
                x, u_up = pycott.cott_utils.get_var_vs_x("u_1", t[i+1], tec_file)      # u1 = u2 = u
                x, rho1_up = pycott.cott_utils.get_var_vs_x("`r_1", t[i+1], tec_file) 
                x, alpha1_up = pycott.cott_utils.get_var_vs_x("`a_1", t[i+1], tec_file)
                x, rho2_up = pycott.cott_utils.get_var_vs_x("`r_2", t[i+1], tec_file)
                x, alpha2_up = pycott.cott_utils.get_var_vs_x("`a_2", t[i+1], tec_file)
                x, Saru_up   = pycott.cott_utils.get_var_vs_x("`S`a`ru", t[i+1], tec_file)
                x, p_up      = pycott.cott_utils.get_var_vs_x("p_mix", t[i+1], tec_file)
                x, gasfrac_up = pycott.cott_utils.get_var_vs_x("gasfrac", t[i+1], tec_file)
                x, mu1_up   = pycott.cott_utils.get_var_vs_x("`m_1", t[i+1], tec_file)
                x, mu2_up   = pycott.cott_utils.get_var_vs_x("`m_2", t[i+1], tec_file)
                x, mu_int_up =  pycott.cott_utils.get_var_vs_x("`m_mix", t[i+1], tec_file)
                x, T1_up     =  pycott.cott_utils.get_var_vs_x("Temp_1", t[i+1], tec_file)
                x, T2_up     =  pycott.cott_utils.get_var_vs_x("Temp_2", t[i+1], tec_file)
                # this is hopefully heat from the wall # no heat in this case (noneq vp)
                #x, Q_w_up = pycott.cott_utils.get_var_vs_x("heat", t[i+1], tec_file)
                # this is hopefully temperature in the wall
                #x, T_w_up  = pycott.cott_utils.get_var_vs_x("T_w", t[i+1], tec_file)  

                if friction:
                    x, fric_up   = pycott.cott_utils.get_var_vs_t("fric", t[i+1], tec_file)
            
                rho_up = rho1_up*alpha1_up + rho2_up*alpha2_up
                rhos_tot_up = np.sum(s_up*rho_up)
                
                # approximate the entropy production
                dt = t[i+1]-t[i]
                drhos_tot_dt = (rhos_tot_up-rhos_tot)/dt # fwd euler in time
                #dt = t[i]-t[i-1]
                #drhos_tot_dt = (rhos_tot-rhos_tot_dwn)/dt # backw euler in time
                # dt = t[i+1]-t[i-1]
                # drhos_tot_dt = (rhos_tot_up-rhos_tot_dwn)/dt # central diff in time
                tot_dSirrdt_num[i] = drhos_tot_dt#*dx #? *dx?
                
                # K_g = noneq_relax_const_K*(mu1-mu2) #hope 1=liq, 2=vap
                # H_g = noneq_relax_const_H*(T1-T2)
                # entr_prod =  H_g*((1/T2)-(1/T1)) + K_g*((mu1-mu_int)/T1 - (mu2-mu_int)/T2)
                K_g = noneq_relax_const_K*(mu2-mu1) #maybe 1=vap, 2=liq?
                H_g = noneq_relax_const_H*(T2-T1)
                entr_prod =  H_g*((1/T1)-(1/T2)) + K_g*((mu2-mu_int)/T2 - (mu1-mu_int)/T1)
                if friction:
                    entr_prod += u*fric*(alpha1*rho1/(rho*T1) + alpha2*rho2/(rho*T2))
                entr_prod_tot = np.sum(entr_prod)
                tot_dSirrdt_calc[i] = entr_prod_tot#*dx
    return t[min_timestep:max_timestep], tot_dSirrdt_num, tot_dSirrdt_calc

def get_cons_vars_cott(tec_file,x, noneq_relax_const_K=0, noneq_relax_const_H=0, friction=False):
    x_from_tec = tec_file.ts[0]['x']
    print("x", x_from_tec)
    dx = x_from_tec[4]-x_from_tec[3]
    
    icell = (np.abs(x_from_tec - x)).argmin()
    
    actual_x = x_from_tec[icell]
    print("Wanted position, x:", x)
    print("Closest position found in tec file, x_tec:", actual_x)
    
    timesteps = tec.nzones #number of time steps in tec-file
    dSirrdt_num = np.zeros(timesteps)
    dSirrdt_calc = np.zeros(timesteps)

    # plot to double check that u are in a rarefaction wave
    u1 = tec.ts[timesteps-1]['u_1']
    # mu = tec.ts[timesteps-1]["`m_1"]
    ##p = tec.ts[timesteps-1]['p']
    plt.plot(x_from_tec,u1)
    plt.xlabel("x")
    plt.ylabel(r"u")
    plt.show()

    # this is just to get the times in the simulation:
    t, hhh     =  pycott.cott_utils.get_var_vs_t("Temp_2", 1, tec_file)
    if (timesteps > 2):
        for i in range(1,timesteps-1):
            # values at t_{i-1}
            x, s_dwn = pycott.cott_utils.get_var_vs_x("entrmix", t[i-1], tec_file)  # = s [J/kg]
            x, u_dwn = pycott.cott_utils.get_var_vs_x("u_1", t[i-1], tec_file)      # u1 = u2 = u
            x, rho1_dwn = pycott.cott_utils.get_var_vs_x("`r_1", t[i-1], tec_file) 
            x, alpha1_dwn = pycott.cott_utils.get_var_vs_x("`a_1", t[i-1], tec_file)
            x, rho2_dwn = pycott.cott_utils.get_var_vs_x("`r_2", t[i-1], tec_file)
            x, alpha2_dwn = pycott.cott_utils.get_var_vs_x("`a_2", t[i-1], tec_file)
            x, p_dwn      = pycott.cott_utils.get_var_vs_x("p_mix", t[i-1], tec_file)
            x, gasfrac_dwn = pycott.cott_utils.get_var_vs_x("gasfrac", t[i-1], tec_file)
            x, mu1_dwn   = pycott.cott_utils.get_var_vs_x("`m_1", t[i-1], tec_file)
            x, mu2_dwn   = pycott.cott_utils.get_var_vs_x("`m_2", t[i-1], tec_file)
            x, mu_int_dwn =  pycott.cott_utils.get_var_vs_x("`m_mix", t[i-1], tec_file)
            x, T1_dwn     =  pycott.cott_utils.get_var_vs_x("Temp_1", t[i-1], tec_file)
            x, T2_dwn     =  pycott.cott_utils.get_var_vs_x("Temp_2", t[i-1], tec_file)
            # this is hopefully heat from the wall # no heat in this case (noneq vp)
            #x, Q_w = pycott.cott_utils.get_var_vs_x("heat", t[i], tec_file)
            # this is hopefully temperature in the wall
            #x, T_w  = pycott.cott_utils.get_var_vs_x("T_w", t[i], tec_file)  
            x, E1_dwn     =  pycott.cott_utils.get_var_vs_x("`a`re_t_1", t[i-1], tec_file) # tot energy for phase 1?
            x, E2_dwn     =  pycott.cott_utils.get_var_vs_x("`a`re_t_2", t[i-1], tec_file) # --"--   for phase 2?
            
            if friction:
                t, fric_dwn   = pycott.cott_utils.get_var_vs_t("fric", t[i-1], tec_file)
            rho_dwn = rho1_dwn*alpha1_dwn + rho2_dwn*alpha2_dwn
                
            # values at t_i
            x, s = pycott.cott_utils.get_var_vs_x("entrmix", t[i], tec_file)  # = s [J/kg]
            x, u = pycott.cott_utils.get_var_vs_x("u_1", t[i], tec_file)      # u1 = u2 = u
            x, rho1 = pycott.cott_utils.get_var_vs_x("`r_1", t[i], tec_file) 
            x, alpha1 = pycott.cott_utils.get_var_vs_x("`a_1", t[i], tec_file)
            x, rho2 = pycott.cott_utils.get_var_vs_x("`r_2", t[i], tec_file)
            x, alpha2 = pycott.cott_utils.get_var_vs_x("`a_2", t[i], tec_file)
            x, p      = pycott.cott_utils.get_var_vs_x("p_mix", t[i], tec_file)
            x, gasfrac = pycott.cott_utils.get_var_vs_x("gasfrac", t[i], tec_file)
            x, mu1   = pycott.cott_utils.get_var_vs_x("`m_1", t[i], tec_file)
            x, mu2   = pycott.cott_utils.get_var_vs_x("`m_2", t[i], tec_file)
            x, mu_int =  pycott.cott_utils.get_var_vs_x("`m_mix", t[i], tec_file)
            x, T1     =  pycott.cott_utils.get_var_vs_x("Temp_1", t[i], tec_file)
            x, T2     =  pycott.cott_utils.get_var_vs_x("Temp_2", t[i], tec_file)
            # this is hopefully heat from the wall # no heat in this case (noneq vp)
            #x, Q_w = pycott.cott_utils.get_var_vs_x("heat", t[i], tec_file)
            # this is hopefully temperature in the wall
            #x, T_w  = pycott.cott_utils.get_var_vs_x("T_w", t[i], tec_file)
            x, E1     =  pycott.cott_utils.get_var_vs_x("`a`re_t_1", t[i], tec_file) # tot energy for phase 1?
            x, E2     =  pycott.cott_utils.get_var_vs_x("`a`re_t_2", t[i], tec_file) # --"--   for phase 2?

            if friction:
                t, fric   = pycott.cott_utils.get_var_vs_t("fric", t[i], tec_file) 
            rho = rho1*alpha1 + rho2*alpha2

            # values at t_{i+1}
            x, s_up = pycott.cott_utils.get_var_vs_x("entrmix", t[i+1], tec_file)  # = s [J/kg]
            x, u_up = pycott.cott_utils.get_var_vs_x("u_1", t[i+1], tec_file)      # u1 = u2 = u
            x, rho1_up = pycott.cott_utils.get_var_vs_x("`r_1", t[i+1], tec_file) 
            x, alpha1_up = pycott.cott_utils.get_var_vs_x("`a_1", t[i+1], tec_file)
            x, rho2_up = pycott.cott_utils.get_var_vs_x("`r_2", t[i+1], tec_file)
            x, alpha2_up = pycott.cott_utils.get_var_vs_x("`a_2", t[i+1], tec_file)
            x, p_up      = pycott.cott_utils.get_var_vs_x("p_mix", t[i+1], tec_file)
            x, gasfrac_up = pycott.cott_utils.get_var_vs_x("gasfrac", t[i+1], tec_file)
            x, mu1_up   = pycott.cott_utils.get_var_vs_x("`m_1", t[i+1], tec_file)
            x, mu2_up   = pycott.cott_utils.get_var_vs_x("`m_2", t[i+1], tec_file)
            x, mu_int_up =  pycott.cott_utils.get_var_vs_x("`m_mix", t[i+1], tec_file)
            x, T1_up     =  pycott.cott_utils.get_var_vs_x("Temp_1", t[i+1], tec_file)
            x, T2_up     =  pycott.cott_utils.get_var_vs_x("Temp_2", t[i+1], tec_file)
            # this is hopefully heat from the wall # no heat in this case (noneq vp)
            #x, Q_w_up = pycott.cott_utils.get_var_vs_x("heat", t[i+1], tec_file)
            # this is hopefully temperature in the wall
            #x, T_w_up  = pycott.cott_utils.get_var_vs_x("T_w", t[i+1], tec_file)
            x, E1_up     =  pycott.cott_utils.get_var_vs_x("`a`re_t_1", t[i+1], tec_file) # tot energy for phase 1?
            x, E2_up     =  pycott.cott_utils.get_var_vs_x("`a`re_t_2", t[i+1], tec_file) # --"--   for phase 2?
            
            if friction:
                x, fric_up   = pycott.cott_utils.get_var_vs_t("fric", t[i+1], tec_file)
            rho_up = rho1_up*alpha1_up + rho2_up*alpha2_up

            dt = t[i+1]-t[i]
            rhous_halfup = dummyflux_LxW(s, rho, u, p, dx, dt, icell)#0.5*(rho[icell+1]*u[icell+1]*s[icell+1]+rho[icell]*u[icell]*s[icell])
            #try simple#MUSCL(conserved_vars, conserved_vars_up)
            rhous_halfdwn = dummyflux_LxW(s, rho, u, p, dx, dt, icell-1)#0.5*(rho[icell]*u[icell]*s[icell]+rho[icell-1]*u[icell-1]*s[icell-1])
            #try simple#MUSCL(conserved_vars_dwn, conserved_vars)
            MUSCL_vars = np.zeros(3)
            MUSCL_vars[0] = rho[icell]
            MUSCL_vars[1] = u[icell]#rho[icell]*u[icell]
            MUSCL_vars[2] = s[icell]#E1[icell]+E2[icell]

            MUSCL_vars_fwd = np.zeros(3)
            MUSCL_vars_fwd[0] = rho[icell+1]
            MUSCL_vars_fwd[1] = u[icell+1]
            MUSCL_vars_fwd[2] = s[icell+1]

            MUSCL_vars_fwdd = np.zeros(3)
            MUSCL_vars_fwdd[0] = rho[icell+2]
            MUSCL_vars_fwdd[1] = u[icell+2]
            MUSCL_vars_fwdd[2] = s[icell+2]

            MUSCL_vars_bcw = np.zeros(3)
            MUSCL_vars_bcw[0] = rho[icell-1]
            MUSCL_vars_bcw[1] = u[icell-1]
            MUSCL_vars_bcw[2] = s[icell-1]

            MUSCL_vars_bcwdd = np.zeros(3)
            MUSCL_vars_bcwdd[0] = rho[icell-2]
            MUSCL_vars_bcwdd[1] = u[icell-2]
            MUSCL_vars_bcwdd[2] = s[icell-2]
            #rhous_halfdwn, rhous_halfup = MUSCL(MUSCL_vars_bcwdd, MUSCL_vars_bcw, MUSCL_vars, MUSCL_vars_fwd, MUSCL_vars_fwdd, dx, dt)
            
            # approximate the entropy production from the simulation
            drhosdt = (rho_up[icell]*s_up[icell]-rho[icell]*s[icell])/dt #fwd Euler in time
            drhousdx= (rhous_halfup-rhous_halfdwn)/dx#dummy_FORCE(s,rho,u,p,dx,dt,icell) #(rhous_half_up-rhous_half_dwn)/dx
            dSirrdt_num[i] = drhosdt + drhousdx
            # calculate theoretical entropy production
            K_g = noneq_relax_const_K*(mu1[icell]-mu2[icell]) #hope 1=liq, 2=vap
            H_g = noneq_relax_const_H*(T1[icell]-T2[icell])
            entr_prod =  H_g*((1/T2[icell])-(1/T1[icell])) + K_g*((mu1[icell]-mu_int[icell])/T1[icell] - (mu2[icell]-mu_int[icell])/T2[icell])
            if friction:
                entr_prod += u[icell]*fric[icell]*(alpha1[icell]*rho1[icell]/(rho[icell]*T1[icell]) + alpha2[icell]*rho2[icell]/(rho[icell]*T2[icell]))
            dSirrdt_calc[i] = entr_prod
    return t[1:-2], dSirrdt_num[1:-2], dSirrdt_calc[1:-2], actual_x

def MUSCL(rhous_bcwdd, rhous_bcw, rhous, rhous_fwd, rhous_fwdd, dx, dt):
    print("rhous, rhousbcw:",rhous, rhous_bcw)
    rhous_L_halfup = rhous + 0.5*minmod((rhous-rhous_bcw)/(rhous_fwd-rhous))*(rhous_fwd-rhous)
    rhous_R_halfup = rhous_fwd + 0.5*minmod((rhous_fwd-rhous)/(rhous_fwdd-rhous_fwd))*(rhous_fwdd-rhous_fwd)
    rhous_L_halfdwn = rhous_bcw + 0.5*minmod((rhous_bcw-rhous_bcwdd)/(rhous-rhous_bcw))*(rhous-rhous_bcw)
    rhous_R_halfdwn = rhous-0.5*minmod((rhous-rhous_bcw)/(rhous_fwd-rhous))*(rhous_fwd-rhous)

    # Kurganov and Tadmor central scheme
    # c_L_halfup = getSpeedOfSound(rhous_L_halfup)
    # c_R_halfup = getSpeedOfSound(rhous_R_halfup)
    # specrad_L_halfup = np.max(np.abs(rhous_L_halfup[1]+c_L_halfup[1]), np.abs(rhous_L_halfup-c_L_halfup))
    # specrad_R_halfup = np.max(np.abs(rhous_R_halfup[1]+c_R_halfup[1]), np.abs(rhous_R_halfup-c_R_halfup))
    # a_halfup = np.max(specrad_L_halfup, specrad_R_halfup)
    a_halfup = dx/dt
    flux_halfup  = 0.5*((rhous_R_halfup[0]*rhous_R_halfup[1]*rhous_R_halfup[2]+
                         rhous_L_halfup[0]*rhous_L_halfup[1]*rhous_L_halfup[1])\
                         - a_halfup*(rhous_R_halfup[2]-rhous_L_halfup[2]))

    # c_L_halfdwn = getSpeedOfSound(rhous_L_halfdwn)
    # c_R_halfdwn = getSpeedOfSound(rhous_R_halfdwn)
    # specrad_L_halfdwn = np.max(np.abs(rhous_L_halfdwn[1]+c_L_halfdwn[1]), np.abs(rhous_L_halfdwn-c_L_halfdwm))
    # specrad_R_halfdwn = np.max(np.abs(rhous_R_halfdwn[1]+c_R_halfdwn[1]), np.abs(rhous_R_halfup-c_R_halfdwn))
    # a_halfdwn = np.max(specrad_L_halfdwn, specrad_R_halfdwn)
    a_halfdwn = dx/dt
    flux_halfdwn = 0.5*((rhous_R_halfdwn[0]*rhous_R_halfdwn[1]*rhous_R_halfdwn[2]+
                        rhous_L_halfdwn[0]*rhous_L_halfdwn[1]*rhous_L_halfdwn[1])\
                        - a_halfdwn*(rhous_R_halfdwn[2]-rhous_L_halfdwn[2]))
    return flux_halfdwn, flux_halfup

def getSpeedOfSound(rhous):
    print("idk what the speed of sound is?")
    return

def minmod(r):
    slope = np.zeros(3)
    if (np.isnan(r[0]) or np.isinf(r[0])): #indicates it is nan
        print("r in minmod:", r)
    else:
        slope[0] = np.maximum(0,np.minimum(1,r[0]))
    if (np.isnan(r[1]) or np.isinf(r[1])): #indicates it is nan
        print("r in minmod:", r)
    else:
        slope[1] = np.maximum(0,np.minimum(1,r[1]))

    if (np.isnan(r[2]) or np.isinf(r[2])): #indicates it is nan
        print("r in minmod:", r)
    else:
        slope[2] = np.maximum(0.0,np.minimum(1.0,r[2]))
    return slope

def dummy_FORCE(s, rho, u, p, dx, dt, j):
    F_iplushalf = 0.5*(dummyflux_LxF(s, rho, u, dx, dt, j)
                       +dummyflux_LxW(s, rho, u, p, dx, dt, j))
    F_iminhalf = 0.5*(dummyflux_LxF(s, rho, u, dx, dt, j-1)
                       +dummyflux_LxW(s, rho, u, p, dx, dt, j-1))
    differential = (F_iplushalf-F_iminhalf)/dx
    return differential

def dummyflux_LxF(s, rho, u, dx, dt, j):
    F_iplushalf = 0.5*(rho[j]*u[j]*s[j]+rho[j+1]*u[j+1]*s[j+1])\
        -0.5*(dx/dt)*(rho[j+1]*s[j+1]-rho[j]*s[j])
    return F_iplushalf

def dummyflux_LxW(s, rho, u, p, dx, dt, j):
    rhos_iplushalf = 0.5*(rho[j]*s[j]+rho[j+1]*s[j+1])\
        -0.5*(dt/dx)*(rho[j+1]*u[j+1]*s[j+1]-rho[j]*u[j]*s[j])
    rho_iplushalf = 0.5*(rho[j]+rho[j+1])\
        -0.5*(dt/dx)*(rho[j+1]*u[j+1]-rho[j]*u[j])
    rhou_iplushalf = 0.5*(rho[j]*u[j]+rho[j+1]*u[j+1])\
        -0.5*(dt/dx)*((rho[j+1]*(u[j+1])**2+p[j+1])-(rho[j]*(u[j])**2+p[j]))
    u_iplushalf = rhou_iplushalf/rho_iplushalf
    
    F_iplushalf = u_iplushalf*rhos_iplushalf
    return F_iplushalf


def plot_specific_entropy_vs_t(tec):
    timesteps = tec.nzones #number of time steps in tec-file
    # hacky way of getting the list of times
    t, hhh     =  pycott.cott_utils.get_var_vs_t("Temp_2", 1, tec)
    stot = np.zeros(timesteps)
    if (timesteps > 2):
        for i in range(0,timesteps):
            x, s     =  pycott.cott_utils.get_var_vs_x("entrmix", t[i], tec)
            stot[i] = np.sum(s)
            
        plt.plot(t, stot)
        plt.title("Total specific entropy in the computational domain vs t")
        plt.ylabel(r"$s$ (Jkg$^{-1}$)")
        plt.xlabel("t (s)")
        plt.show()
        
        dsdt_tot = np.gradient(stot,t)
        plt.plot(t, dsdt_tot)
        plt.title("Total change of specific entropy over time in the computational domain")
        plt.ylabel(r"$\frac{\partial s}{\partial t}$ (Jkg$^{-1}$s$^{-1}$)")
        plt.xlabel("t (s)")
        plt.show()
    else:
        print("Not enough timesteps available to find dsdt")
    return

def plot_specific_entropy_vs_t_CV(tec, x):
    x_from_tec = tec.ts[0]['x']
    icell = (np.abs(x_from_tec - x)).argmin()    
    actual_x = x_from_tec[icell]
    print("Wanted position, x:", x)
    print("Closest position found in tec file, x_tec:", actual_x)
    
    t, s     =  pycott.cott_utils.get_var_vs_t("entrmix", icell, tec)
    t, rho1 = pycott.cott_utils.get_var_vs_t("`r_1", icell, tec) 
    t, alpha1 = pycott.cott_utils.get_var_vs_t("`a_1", icell, tec)
    t, rho2 = pycott.cott_utils.get_var_vs_t("`r_2", icell, tec)
    t, alpha2 = pycott.cott_utils.get_var_vs_t("`a_2", icell, tec)
    rho = rho1*alpha1 + rho2*alpha2
            
    plt.plot(t, s)
    plt.title("Specific entropy in CV w cell center at x= "+ str(actual_x))
    plt.ylabel(r"$s$ (Jkg$^{-1}$)")
    plt.xlabel("t (s)")
    plt.show()

    plt.plot(t, s*rho)
    plt.title(r"$\rho s$ in CV w cell center at x= "+ str(actual_x))
    plt.ylabel(r"$\rho s$ (Jm$^{-3}$)")
    plt.xlabel("t (s)")
    plt.show()

    dsdt = np.gradient(s,t)
    plt.plot(t, dsdt)
    plt.title(r"Total change of $s$ over time in CV w cell center at x= "+ str(actual_x))
    plt.ylabel(r"$\frac{\partial s}{\partial t}$ (Jkg$^{-1}$s$^{-1}$)")
    plt.xlabel("t (s)")
    plt.show()
    
    drhosdt = np.gradient(s*rho,t)
    plt.plot(t, drhosdt)
    plt.title(r"Total change of $\rho s$ over time in CV w cell center at x= "+ str(actual_x))
    plt.ylabel(r"$\frac{\partial\rho s}{\partial t}$ (Jm$^{-3}$s$^{-1}$)")
    plt.xlabel("t (s)")
    plt.show()
    return

#------------------------------Full domain calculation w Gauss curve-------------------------------------
# Simulation gives negative entropy production??

# casedir  = "../Simulations/noneq_gausscurve"
# tec_path = os.path.join(casedir, "cott.tec")
# inp_path = os.path.join(casedir, "user.inp")
# tec =  pycott.cott_utils.load_tec(tec_path)  # tec object for the simulation
# times1, sigma1, sigma1calc = get_cott_sigma_vs_t_at_x(tec, 1, 40, noneq_relax_const_K=0.0,
#                                                                    noneq_relax_const_H=0.0)
# plt.plot(times1, sigma1, label="Simulated: $\mathcal{K}=0, \mathcal{H}=0$", linestyle="--")
# plt.plot(times1, sigma1calc, label="Analytic: $\mathcal{K}=0, \mathcal{H}=0$")
# plt.title(r"$\frac{dSirr}{dt}$ in calculation domain (0m-120m) over time ")
# plt.xlabel("t (s)")
# plt.ylabel(r"$\frac{dSirr}{dt}$ (J K$^{-1}$ m$^{-2}$ s$^{-1}$)")
# plt.legend()
# plt.tight_layout()
# plt.show()

# casedir  = "../Simulations/noneq_gausscurve5000"
# tec_path = os.path.join(casedir, "cott.tec")
# inp_path = os.path.join(casedir, "user.inp")
# tec =  pycott.cott_utils.load_tec(tec_path)  # tec object for the simulation
# times1, sigma1, sigma1calc = get_cott_sigma_vs_t_at_x(tec, 1, 40, noneq_relax_const_K=0.0,
#                                                                    noneq_relax_const_H=0.0)
# plt.plot(times1, sigma1, label="Simulated: $\mathcal{K}=0, \mathcal{H}=0$", linestyle="--")
# plt.plot(times1, sigma1calc, label="Analytic: $\mathcal{K}=0, \mathcal{H}=0$")
# plt.title(r"$\frac{dSirr}{dt}$ in calculation domain (0m-120m) over time ")
# plt.xlabel("t (s)")
# plt.ylabel(r"$\frac{dSirr}{dt}$ (J K$^{-1}$ m$^{-2}$ s$^{-1}$)")
# plt.legend()
# plt.tight_layout()
# plt.show()

# casedir  = "../Simulations/noneq_gausscurve"
# tec_path = os.path.join(casedir, "cott.tec")
# inp_path = os.path.join(casedir, "user.inp")
# tec =  pycott.cott_utils.load_tec(tec_path)  # tec object for the simulation
# plot_specific_entropy_vs_t(tec)

# plot_specific_entropy_vs_t_CV(tec, 20)


#------------------------------Full domain calculation w shock-------------------------------------
# analytic doesn't take shock into account :C? does not give right answer

# casedir  = "../Simulations/entropy_production_vp1000"
# tec_path = os.path.join(casedir, "cott.tec")
# inp_path = os.path.join(casedir, "user.inp")
# tec =  pycott.cott_utils.load_tec(tec_path)  # tec object for the simulation
# times, sigma, sigmacalc = get_cott_sigma_vs_t_at_x(tec, 1, 100)
# plt.plot(times, sigma, label="Simulated: $\mathcal{K}=0, \mathcal{H}=0$", linestyle="--")
# plt.plot(times, sigmacalc, label="Analytic: $\mathcal{K}=0, \mathcal{H}=0$")
# plt.title(r"$\frac{dSirr}{dt}$ in calculation domain (0m-120m) over time ")
# plt.xlabel("t (s)")
# plt.ylabel(r"$\frac{dSirr}{dt}$ (J K$^{-1}$ m$^{-2}$ s$^{-1}$)")
# plt.legend()
# plt.tight_layout()
# plt.show()

# casedir  = "../Simulations/entropy_production_vp5000"
# tec_path = os.path.join(casedir, "cott.tec")
# inp_path = os.path.join(casedir, "user.inp")
# tec =  pycott.cott_utils.load_tec(tec_path)  # tec object for the simulation
# times, sigma, sigmacalc  = get_cott_sigma_vs_t_at_x(tec, 1, 100)
# plt.plot(times, sigma, label="Simulated: $\mathcal{K}=0, \mathcal{H}=0$", linestyle="--")
# plt.plot(times, sigmacalc, label="Analytic: $\mathcal{K}=0, \mathcal{H}=0$")
# plt.title(r"$\frac{dSirr}{dt}$ in calculation domain (0m-120m) over time ")
# plt.xlabel("t (s)")
# plt.ylabel(r"$\frac{dSirr}{dt}$ (J K$^{-1}$ m$^{-2}$ s$^{-1}$)")
# plt.legend()
# plt.tight_layout()
# plt.show()

# ---------------------------- Control volume calculation w shock-----------------------------------
# gives positive then negative entropy production as rarefaction wave passes through gridcell?
# should get lower for finer grid if this is caused by discretization error on my side

#------------------No source----------------------
# No source = theoretically no entropy production
casedir  = "../Simulations/entropy_production_vp1000"
tec_path = os.path.join(casedir, "cott.tec")
inp_path = os.path.join(casedir, "user.inp")
tec =  pycott.cott_utils.load_tec(tec_path)  # tec object for the simulation
times, sigma, sigmacalc, actual_x = get_cons_vars_cott(tec, 70, noneq_relax_const_K=0.0,
                                                                     noneq_relax_const_H=0.0)
#plot_specific_entropy_vs_t_CV(tec, 70)

# # No source = theoretically no entropy production, finer grid
# casedir  = "../Simulations/entropy_production_vp5000"
# tec_path = os.path.join(casedir, "cott.tec")
# inp_path = os.path.join(casedir, "user.inp")
# tec =  pycott.cott_utils.load_tec(tec_path)  # tec object for the simulation
# times1, sigma1, sigmacalc1, actual_x1 = get_cons_vars_cott(tec, 70, noneq_relax_const_K=0.0,
#                                                                      noneq_relax_const_H=0.0)

# plt.plot(times, sigma, label="Simulated: $\mathcal{K}=0, \mathcal{H}=0, N=1000$", linestyle="--")
# plt.plot(times, sigmacalc, label="Analytic: $\mathcal{K}=0, \mathcal{H}=0, N=1000$")
# plt.plot(times1, sigma1, label="Simulated: $\mathcal{K}=0, \mathcal{H}=0, N=5000$", linestyle="--")
# plt.plot(times1, sigmacalc1, label="Analytic: $\mathcal{K}=0, \mathcal{H}=0, N=5000$")
# plt.title(r"$\frac{dSirr}{dt}$ at x = "+  str(actual_x))
# plt.xlabel("t (s)")
# plt.ylabel(r"$\frac{dSirr}{dt}$ (J K$^{-1}$ m$^{-2}$ s$^{-1}$)")
# plt.legend()
# plt.tight_layout()
# plt.show()

#plot_specific_entropy_vs_t_CV(tec, 70)

#------------------Source KE-4 HE4-------------------
casedir  = "../Simulations/noneq_constSource"
tec_path = os.path.join(casedir, "cott.tec")
inp_path = os.path.join(casedir, "user.inp")
tec =  pycott.cott_utils.load_tec(tec_path)  # tec object for the simulation
timesK4H4, sigmaK4H4, sigmacalcK4H4, actual_xK4H4 = get_cons_vars_cott(tec, 70, noneq_relax_const_K=0.0001,
                                                                     noneq_relax_const_H=10000)

# plot_specific_entropy_vs_t_CV(tec, 70)
plt.plot(times, sigma, label="Simulated: $\mathcal{K}=0, \mathcal{H}=0$", linestyle="--")
plt.plot(timesK4H4, sigmaK4H4, label="Simulated: $\mathcal{K}=10^{-4}, \mathcal{H}=10^{4}$", linestyle="--")
plt.plot(timesK4H4, sigmacalcK4H4, label="Analytic: $\mathcal{K}=10^{-4}, \mathcal{H}=10^{4}$")
plt.title(r"$\frac{dSirr}{dt}$ at x = "+  str(actual_x))
plt.xlabel("t (s)")
plt.ylabel(r"$\frac{dSirr}{dt}$ (J K$^{-1}$ m$^{-2}$ s$^{-1}$)")
plt.legend()
plt.tight_layout()
plt.show()

# casedir  = "../Simulations/noneq_constSourceK01HE4"
# tec_path = os.path.join(casedir, "cott.tec")
# inp_path = os.path.join(casedir, "user.inp")
# tec =  pycott.cott_utils.load_tec(tec_path)  # tec object for the simulation
# # times2, sigma2, sigma2calc, actual_x_a2, actual_x_b2 = get_cott_sigma_vs_t_at_x(tec, 0.1, 120,
# #                                                                  noneq_relax_const_K=0.01,
# #                                                                  noneq_relax_const_H=10000)
# times2, sigma2, sigma2calc, actual_x = get_cons_vars_cott(tec, 70,
#                                                                   noneq_relax_const_K=0.01,
#                                                                   noneq_relax_const_H=10000)
                                                                                
# casedir  = "../Simulations/noneq_constSourceKE-4HE8"
# tec_path = os.path.join(casedir, "cott.tec")
# inp_path = os.path.join(casedir, "user.inp")
# tec =  pycott.cott_utils.load_tec(tec_path)  # tec object for the simulation
# # times3, sigma3, sigma3calc, actual_x_a3, actual_x_b3 = get_cott_sigma_vs_t_at_x(tec, 0.1, 120,
# #                                                                  noneq_relax_const_K=0.0001,
# #                                                                  noneq_relax_const_H=1000000)
# times3, sigma3, sigma3calc, actual_x_a3 = get_cons_vars_cott(tec, 70,
#                                                                  noneq_relax_const_K=0.0001,
#                                                                  noneq_relax_const_H=1000000)





# plt.plot(times, sigma, label="Simulated: $\mathcal{K}=0, \mathcal{H}=0$", linestyle="--")
# plt.plot(times2, sigma2, label="Simulated: $\mathcal{K}=10^{-2}, \mathcal{H}=10^{4}$", linestyle="--")
# plt.plot(times2, sigma2calc, label="Analytic: $\mathcal{K}=10^{-2}, \mathcal{H}=10^{4}$")
# plt.title(r"$\frac{dSirr}{dt}$ at x = "+  str(actual_x))
# plt.xlabel("t (s)")
# plt.ylabel(r"$\frac{dSirr}{dt}$ (J K$^{-1}$ m$^{-2}$ s$^{-1}$)")
# plt.legend()
# plt.tight_layout()
# plt.show()

# plt.plot(times, sigma, label="Simulated: $\mathcal{K}=0, \mathcal{H}=0$", linestyle="--")
# plt.plot(times3, sigma3, label="Simulated: $\mathcal{K}=10^{-4}, \mathcal{H}=10^{6}$", linestyle="--")
# plt.plot(times3, sigma3calc, label="Analytic: $\mathcal{K}=10^{-4}, \mathcal{H}=10^{6}$")
# plt.title(r"$\frac{dSirr}{dt}$ at x = "+  str(actual_x))
# plt.xlabel("t (s)")
# plt.ylabel(r"$\frac{dSirr}{dt}$ (J K$^{-1}$ m$^{-2}$ s$^{-1}$)")
# plt.legend()
# plt.tight_layout()
# plt.show()

