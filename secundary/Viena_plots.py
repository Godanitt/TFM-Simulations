import uproot
import os
import matplotlib.pyplot as plt
import numpy as np

#Carpeta donde están guardados todos los outputs de root (directorio de trabajo de python + rootArchives)
base_folder = ''


ar_cf4_experimental_concs = np.array([5, 10, 67, 100])
ar_cf4_experimental_phe_e = np.array([0.575117370892018, 0.548826291079812, 0.300156494522691, 0.0926277372262773])
E_fields_ar_cf4 = {5:65, 10:78, 67:88, 100:95}

he_cf4_experimental_concs = np.array([20, 40, 100])
he_cf4_experimental_phe_e = np.array([0.0668231611893583, 0.0701095461658841, 0.0926277372262773])
E_fields_he_cf4 = {20: 60, 40:75, 100:95}

fit_parameters_wo_uncertainties = { 'f_CF3':0.0963, 'f_ar_exc*P2A1':1, 'quenching_ratio_ar_ex':0.038,
                   'tau_1_K_cool':0.026,  'K_scint_tau_2':0.070, 'f_CF4':0.291,
                    'quenching_ratio_ar_ion':6.238, 'P_CD':0.671, 'tau3_cont':5,
}




def open_root(P, gas_a, conc_a, conc_b, E_field, gap, npe, gas_b='cf4', **kwargs):

    if conc_b == 100:
        root_path = base_folder + f'{gas_b}0.0ar_{E_field:.1f}kVcm_{P:d}bar_{gap:.2f}cm_{npe:d}npe.root'

    else:
        root_path = base_folder + f'{gas_a}{conc_b:.1f}{gas_b}_{E_field:.1f}kVcm_{P:d}bar_{gap:.2f}cm_{npe:d}npe.root'

    root_file = uproot.open(root_path)

    hLevels = root_file['hLevels'].values()

    return root_file, hLevels




def ph_e(P, gas_a, conc_a, conc_b, E_field, gap, npe, n_start_ar=0, n_start_he=0, gas_b='cf4',**kwargs):

    root_file, hLevels = open_root(P, gas_a, conc_a, conc_b, E_field, gap, npe, gas_b, **kwargs)


    n_e = sum(root_file['dataPerPrimaryElectron;1']['nElectrons'].array())

    n_exc_ar_vis, n_exc_he_vis, n_ion_ar_uv, n_ion_he_uv = 0, 0, 0, 0
    n_exc_vis, n_ion_uv, n_exc_dir_vis, n_ion_dir_uv = 0, 0, 0, 0

    if conc_b == 100:
        
            n_exc_dir_vis = sum(hLevels[19 + 16:])#Teño que empezar no 16º nivel para que me axuste o punto de CF4 puro

            n_ion_dir_uv = hLevels[1] + hLevels[6] #Canal directo es el que produce CF3+ de algún modo

            

    else:
        
        if gas_a == 'ar':
            index_start_ar = 3 + n_start_ar #seleccionamos desde que estado empezamos a contar
            n_exc_ar_vis = sum(hLevels[index_start_ar:47])
            n_exc_vis = n_exc_ar_vis
            n_exc_dir_vis = sum(hLevels[66 + 16:])

            n_ion_ar_uv = 0 #non encontro ionizacións dobles e triples, a sección eficaz debe ser baixa
            n_ion_uv = n_ion_ar_uv
            n_ion_dir_uv = hLevels[48] + hLevels[53] #Canal directo es el que produce CF3+ de algún modo

        elif gas_a == 'he':
            index_start_he = 3 + n_start_he
            n_exc_he_vis = sum(hLevels[index_start_he:52])
            n_exc_vis = n_exc_he_vis
            n_exc_dir_vis = sum(hLevels[71 + 16:])

            n_ion_he_uv = 0 #non encontro ionizacións dobles e triples, a sección eficaz debe ser baixa
            n_ion_uv = n_ion_he_uv
            n_ion_dir_uv = hLevels[53] + hLevels[58] #Canal directo es el que produce CF3+ de algún modo

        else:
            print(f'El gas A introducido: {gas_a} no es válido')

    conc_b *= 1e-2
    conc_a *= 1e-2

    ph_e_vis = n_exc_dir_vis / n_e * kwargs['f_CF3'] + n_exc_vis / n_e * kwargs['f_ar_exc*P2A1'] * (1 / (
            1 + conc_a / conc_b * kwargs['quenching_ratio_ar_ex']) )
    
    ph_e_uv = (n_ion_dir_uv / n_e * kwargs['f_CF4']+ n_ion_uv / n_e * ( P * conc_b * kwargs['quenching_ratio_ar_ion'] / ( 
            1 / kwargs['tau3_cont']  + P * conc_b * kwargs['quenching_ratio_ar_ion']) ) * kwargs['P_CD']) * (
            P * conc_b / (P * conc_b + kwargs['tau_1_K_cool']
            )  * 1 / (1 + kwargs['K_scint_tau_2'] * P * conc_b))
    
    dic_values = {
        'n_exc_ar_vis':n_exc_ar_vis, 'n_exc_he_vis':n_exc_he_vis, 'n_ion_ar_uv':n_ion_ar_uv, 
        'n_ion_he_uv':n_ion_he_uv, 'n_exc_dir_vis':n_exc_dir_vis, 'n_ion_dir_uv':n_ion_dir_uv,
    }

    return ph_e_vis, ph_e_uv, dic_values




figure_ar, ax_ar = plt.subplots(dpi=200)
ax_ar.plot(ar_cf4_experimental_concs / 100, ar_cf4_experimental_phe_e, 'ks', label='Measured points (LIP group)')
ax_ar.set_title('Argon mixtures')
ax_ar.set_xlim(0, 1.05)
#ax_ar.set_ylim(bottom=0)
ax_ar.set_xticks(np.arange(0, 1.1, 0.1))
ax_ar.set_xlabel('CF4 Fraction')
ax_ar.set_ylabel('ph/e$^-$')


figure_he, ax_he = plt.subplots(dpi=200)
ax_he.plot(he_cf4_experimental_concs / 100, he_cf4_experimental_phe_e, 'ks')
ax_he.set_title('Helium mixtures')
ax_he.set_xlim(0, 1.05)

ax_he.set_xticks(np.arange(0, 1.1, 0.1))
ax_he.set_xlabel('CF4 Fraction')
ax_he.set_ylabel('ph/e$^-$')


# ph_e_ar_plot = []
# for conc in ar_cf4_experimental_concs:
#     E_field = E_fields_ar_cf4[conc]
#     ph_e_vis_cf4, ph_e_uv_cf4, dic_values_cf4 = ph_e(1, 'ar', 1-conc, conc, E_field, 0.05, 10000, n_start_ar=29, **fit_parameters_wo_uncertainties)

#     print(dic_values_cf4)

#     ph_e_ar_plot.append(ph_e_vis_cf4)

# ax_ar.plot(ar_cf4_experimental_concs / 100, ph_e_ar_plot, '.', label='n=29')



# ph_e_he_plot = []
# for conc in he_cf4_experimental_concs:
#     E_field = E_fields_he_cf4[conc]
#     ph_e_vis_cf4, ph_e_uv_cf4, dic_values_cf4 = ph_e(1, 'he', 100-conc, conc, E_field, 0.05, 10000, **fit_parameters_wo_uncertainties)

#     print(dic_values_cf4)

#     ph_e_he_plot.append(ph_e_vis_cf4)

# ax_he.plot(he_cf4_experimental_concs/100, ph_e_he_plot, '.', label='n=0')




dic_states_ar = {18:'EXC 3D4    J=3                    ELEVEL= 14.013 ',
              22:'EXC 3D1!   J=3                    ELEVEL= 14.099 ',
              26:'EXC 3S1!!! J=3                    ELEVEL= 14.236 ',
              29:' EXC 3S1!   J=1 RESONANT           ELEVEL= 14.304'
              }

dic_states_ar = {18:'14.013 eV',
              22:'14.099 eV',
              26:'14.236 eV',
              29:'14.304 eV'
              }

n_list = [18,22,26,29]

for n in n_list:
    ph_e_ar_plot = []
    for conc in ar_cf4_experimental_concs:
        E_field = E_fields_ar_cf4[conc]
        ph_e_vis_cf4, ph_e_uv_cf4, dic_values_cf4 = ph_e(1, 'ar', 1-conc, conc, E_field, 0.05, 10000, n_start_ar=n, **fit_parameters_wo_uncertainties)

        #print(dic_values_cf4)

        ph_e_ar_plot.append(ph_e_vis_cf4)
    if n in dic_states_ar.keys():
        label = dic_states_ar[n]
    else:
        label=None
    ax_ar.plot(ar_cf4_experimental_concs / 100, ph_e_ar_plot, '.', label=label)
ax_ar.legend()


for n in range(0, 50):
    ph_e_he_plot = []
    for conc in he_cf4_experimental_concs:
        E_field = E_fields_he_cf4[conc]
        ph_e_vis_cf4, ph_e_uv_cf4, dic_values_cf4 = ph_e(1, 'he', 100-conc, conc, E_field, 0.05, 10000, n_start_he=n,**fit_parameters_wo_uncertainties)

        #print(dic_values_cf4)

        ph_e_he_plot.append(ph_e_vis_cf4)
    ax_he.plot(he_cf4_experimental_concs/100, ph_e_he_plot, '.', label=f'n={n}')
#ax_he.legend()
ax_he.set_ylim(bottom=0)

figure_ar.savefig("Secundario.pdf",bbox_inches="tight")