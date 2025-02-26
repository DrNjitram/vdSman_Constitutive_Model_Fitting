import os
import re
from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import namedtuple

from matplotlib.patches import Patch
from scipy.optimize import curve_fit

Measurement = namedtuple("Measurement", "Mixture, Customer, Instrument, Material, Comments, Info, Batch_NB, Order_NB, N3, N4, N5, Frequency, Temperature, InterID, Moisture, pH, TgT, Plateau, Plateau_std, tan_delta, data")
#path = r"D:\Data\Rheology\CCR\202409 DOE pH T M Arcon F\Overview.xlsx"
colors = ["k", "r", "g", "b", "m", "c", "y", "k", "k", "k", "k", "k", "k"]
Cpw = 1.92E+03
Tgw = 139 #K
Cps = 4.25E+02
Tgs = 384 #K
debug = False
All_measurements = dict()
unique_id = []
TgT_eq = lambda x, *params: params[0] + params[1] / (1 + np.pow(x / params[2], params[3]))
linear_eq = lambda x, *params: params[0] + params[1]*x
# pH correction data
# 2nd order poly
base_data = [1.7891, 3.9621, 6.3907]
acid_data = [4.2209, -5754, 6.4553]

def fit_equation(x_vals, y_vals, linear=False):
    if linear:
        params = [0, 1]
        params, covariance = curve_fit(linear_eq, x_vals, np.log10(y_vals), p0=params)
    else:
        params = [1500, -2000, 0.5, 9]
        params, covariance = curve_fit(TgT_eq, x_vals, y_vals, p0=params)
    return params

def TgT(M, T):
    M = float(M)
    if type(T) is str and "C" in T:
        T = T[:-2]
    T = float(T)
    return ((M * Cpw * Tgw + (1 - M) * Cps * Tgs) / (M * Cpw + (1 - M) * Cps))/(T + 273)

def backcalculate_moisture(_TgT, T):
    return (Cps * Tgs - _TgT * (T + 273) * Cps)/(_TgT * (T + 273) * Cpw - _TgT * (T + 273) * Cps - Cpw * Tgw + Cps * Tgs)

def recalculate_pH(pH):
    if pH < 7.0:
        x2, x1, x0 = acid_data
        A = -(10**-pH)
    else:
        x2, x1, x0 = base_data
        A = 10 ** -(14-pH)
    return x2*(A**2)+x1*A+x0

# dirty excel regression
excel_regr = lambda x: -159.13*x**4 - 44.344*x**3+24.863*x**2 + 13.665*x+6.6294

def get_measurements(path, filenames, measurement_width, header_row, headers, default_pH=7.0, re_calculate_pH=False, mixture=None) -> List[Measurement]:
    measurements = []
    for filename in filenames:
        index = 0
        file_path = os.path.join(path, filename)
        print(filename)
        while True:
            try:
                info = pd.read_excel(io=file_path, header=None, skiprows=1, nrows=14, usecols=[(1 if not index else 0)+index*measurement_width], names=["Value"])
                data = pd.read_excel(io=file_path, header=None, names=headers, skiprows=header_row, usecols=list(range(index * measurement_width, (index+1)*measurement_width)))
                measurement = Measurement(*info["Value"], -1, default_pH, -1, -1, -1, -1, data)
                if measurement.InterID not in unique_id:
                    unique_id.append(measurement.InterID)
                    information_string = f"{measurement.Material} {measurement.Info} {measurement.Comments}"
                    for i, bit in enumerate(information_string.split()):
                        if ("M" in bit and "%" in bit) or bit.startswith("M"):
                            measurement = measurement._replace(Moisture = float(re.sub(r"\D", "", bit))/100)
                        elif "pH" in bit: # Old pH stuff
                            if len(bit) == 2:
                                pH = float(information_string.split()[i+1])
                            else:
                                pH = float(bit[2:])

                            if re_calculate_pH:
                                measurement = measurement._replace(pH=recalculate_pH(pH))
                            else:
                                measurement = measurement._replace(pH=float(pH))
                        elif bit.endswith("M") or "." in bit:
                            bit = float(re.sub(r"[^0-9.,]", "", bit))
                            if "HCl" in information_string:
                                pH = -1*bit
                            elif "NaOH" in information_string:
                                pH = bit
                            elif bit != 0.0:
                                if bit == 0.752:
                                    pH = 0.52
                                elif bit in [0.07, 0.24, 0.41, 0.74, 0.3]:
                                    pH = -1*bit
                                else:
                                    pH = bit
                            else:
                                pH = 0.0


                            measurement = measurement._replace(pH=float(excel_regr(pH*measurement.Moisture)))
                    if measurement.N3 != 0:
                        measurement = measurement._replace(pH=measurement.N3)
                    if measurement.Material == "A8" and measurement.pH == 0.0:
                        measurement = measurement._replace(pH=float(excel_regr(0)))
                    if mixture:
                        measurement = measurement._replace(Material=mixture)




                    measurement = measurement._replace(Temperature=float(measurement.Temperature[:-2]))
                    measurement = measurement._replace(Plateau=np.mean(measurement.data["G' [kPa]"][5:10]))
                    measurement = measurement._replace(Plateau_std=np.std(measurement.data["G' [kPa]"][5:10]))
                    measurement = measurement._replace(TgT = TgT(measurement.Moisture, measurement.Temperature))
                    measurement = measurement._replace(tan_delta=np.mean(measurement.data["Tan Delta"][5:10]/measurement.Plateau))

                    measurements.append(measurement)
                    print(f"{len(measurements)-1}: {measurement.Plateau:.2f}, {measurement[:-1]}")
                    if  debug and measurement.Temperature == 150.0 and measurement.Moisture == 0.55:

                        plt.figure()
                        plt.loglog(measurement.data["Strain %"], measurement.data["G' [kPa]"], f"ko", label=measurement.Info)
                        plt.loglog(measurement.data["Strain %"], measurement.data["G'' [kPa]"], f"ks")
                        plt.plot(np.mean(measurement.data["Strain %"][5:10]), measurement.Plateau, "ro")
                        plt.title(f"{len(measurements)}, {measurement.Temperature}, {measurement.Moisture}, {measurement.pH}, {measurement.TgT:.2f},  {measurement.Plateau:.2f}")
                        plt.show()

                index += 1

            except Exception as e:
                if not measurements:
                    print(e)
                    raise e
                break
    return measurements

def graph_measurement(measurement):
    plt.loglog(measurement.data["Strain %"], measurement.data["G' [kPa]"], f"{colors[0]}o", label="G'")
    plt.loglog(measurement.data["Strain %"], measurement.data["G'' [kPa]"], f"{colors[0]}s", label="G\"")
    plt.ylabel("G' [kPa]/G\" [kPa] ")
    plt.xlabel("Strain %")
    plt.legend()
    plt.title(
        f"T{measurement.Temperature}°C, M{measurement.Moisture * 100:.0f}%, pH{measurement.pH:.2f}")
    plt.show()


def ArconF_T():
    header_row = 18

    headers = ['Strain %', "G' [kPa]", "G'' [kPa]", 'G* [kPa]', 'Tan Delta',
               'Temperature [°C]', 'Pressure [bar]', "Viscosity n' [Pa*s]", 'Viscosity n" [Pa*s]',
               'Viscosity n* [Pa*s]', 'Shear Rate [1/sec]', 'Shear Stress [Pa]']

    path = rf"{Drive}:\Data\Rheology\CCR\202409 DOE pH T M Arcon F\Redo"
    filenames = [f"Block_{i}.xls" for i in range(1,6)]
    measurements = get_measurements(path, filenames, len(headers), header_row, headers, 6.42, True, mixture="Arcon F")



    # headers = ['Strain %', "S' [dNm]", "S'' [dNm]", 'S* [dNm]', "G' [kPa]", "G'' [kPa]", 'G* [kPa]', 'Tan Delta',
    #            'Temperature [°C]', 'Pressure [bar]', "Viscosity n' [Pa*s]", 'Viscosity n" [Pa*s]',
    #            'Viscosity n* [Pa*s]', 'Shear Rate [1/sec]', 'I2 abs.', 'I3 abs.', 'I2/I1 [%]', 'I3/I1 [%]',
    #            'Shear Stress [Pa]']
    #
    # path = rf"{Drive}:\Data\Rheology\CCR\202409 DOE pH T M Arcon F"
    # filenames = ["Overview_fixed.xls", "Overview2.xls", "Overview3.xls", "Overview4.xls", ]
    # measurements = get_measurements(path, filenames, len(headers), header_row, headers, 6.42, True, mixture="Arcon F")
    #
    #
    # measurement_width = 12
    # headers = ['Strain %', "G' [kPa]", "G'' [kPa]", 'G* [kPa]', 'Tan Delta',
    #            'Temperature [°C]', 'Pressure [bar]', "Viscosity n' [Pa*s]", 'Viscosity n" [Pa*s]',
    #            'Viscosity n* [Pa*s]', 'Shear Rate [1/sec]', 'Shear Stress [Pa]']
    # measurements2 = get_measurements(path, ["Overview5.xlsx", "Overview6.xlsx"], measurement_width, header_row, headers, 6.42, False)
    #
    # for measurement, pH in zip(measurements2[1:], measured_pH[1:]): # First one is busted
    #     #print(f"{pH}, {measurement.Comments}, {measurement.Info}, {measurement.Temperature}")
    #     measurement = measurement._replace(pH=pH)
    #     measurements.append(measurement)

    pH_InterIDs = [23627, 23630, 23633, 23636, 23639, 23642, 23645, 23648, 23651, 23654, 23657, 23660, 23663, 23666,
                   23669, 23672, 23675, 23678]
    measured_pH = [4.80, 5.01, 5.87, 6.02, 4.78, 5.61, 5.22, 5.17, 5.60, 6.31, 6.46, 6.45, 6.42, 6.85, 8.62, 9.51, 5.03,
                   6.26]
    for i, measurement in enumerate(measurements):
        if measurement.InterID == 23462:
            measurements[i] = False
        if measurement.InterID in pH_InterIDs:
            measurements[i] = measurement._replace(pH=measured_pH[pH_InterIDs.index(measurement.InterID)])



    All_measurements["ArconF"] = measurements


def new_Alpha8():
    measurement_width = 12
    header_row = 18
    headers = ['Strain %', "G' [kPa]", "G'' [kPa]", 'G* [kPa]', 'Tan Delta', 'Temperature [°C]', 'Pressure [bar]', "Viscosity n' [Pa*s]", 'Viscosity n" [Pa*s]', 'Viscosity n* [Pa*s]', 'Shear Rate [1/sec]', 'Shear Stress [Pa]']
    path = rf"{Drive}:\Data\Rheology\CCR\202411 DoE pH T M Alph 8"
    filenames = ["Overview 1.xls", "Overview 2.xls", "Overview 3.xls", "Overview 4.xls"]
    measurements = get_measurements(path, filenames, measurement_width, header_row, headers, default_pH=0, mixture="Alpha 8")
    measurements[7] = False
    measurements[17] = False
    measurements[25] = False
    All_measurements["Alpha8"] = measurements

    measurement_width = 22
    header_row = 18
    headers = ['Strain %', "S' [dNm]", "S'' [dNm]", 'S* [dNm]', "G' [kPa]", "G'' [kPa]", 'G* [kPa]', 'Tan Delta',
               'Temperature [°C]', 'Pressure [bar]', "Viscosity n' [Pa*s]", 'Viscosity n" [Pa*s]',
               'Viscosity n* [Pa*s]', 'Shear Rate [1/sec]', "J' [1/kPa]", "J'' [1/kPa]", "	J* [1/kPa]"
, 'I2 abs.', 'I3 abs.', 'I2/I1 [%]', 'I3/I1 [%]',
               'Shear Stress [Pa]']

    filenames = ["Overview 5.xls"]
    measurements = get_measurements(path, filenames, measurement_width, header_row, headers, 6.42, True,
                                    mixture="Alpha 8")
    All_measurements["Alpha8"] =  All_measurements["Alpha8"]+ measurements


def JMP_export():
    header_string = "ID, Material, T (K), M (%), pH, TgT (-), G0 [kPa], tan(δ)"
    print_measurement = lambda m: print(f"{m.InterID}, {m.Material}, {m.Temperature + 273.15}, {m.Moisture:.2f}, {m.pH:.2f}, {m.TgT:.2f}, {m.Plateau:.2f}, {m.tan_delta:.4f}")

    print(header_string)
    if "ArconF" in All_measurements:
        for measurement in All_measurements["ArconF"]:
            if measurement:
                print_measurement(measurement)
            else:
                print("nan")

    if "Alpha8" in All_measurements:
        for i, measurement in enumerate(All_measurements["Alpha8"]):
            if measurement:
                print_measurement(measurement)
            else:
                print("nan")



def TgT_comparison():
    x_indices = np.linspace(0.38, 0.55, 20)

    x_vals = []
    y_vals = []
    fig, ax = plt.subplots(dpi=600)
    for measurement in All_measurements["Alpha8"]:
        if measurement and 6 < measurement.pH < 7:
            x_vals.append(measurement.TgT)
            y_vals.append(measurement.Plateau)
            plt.errorbar(measurement.TgT, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"ks")

    # params = fit_equation(x_vals, y_vals, True)
    # points = linear_eq(x_indices, *params)
    # plt.plot(x_indices, 10 ** points, "k--", label="fit", zorder=10)
    ax.set_yscale("log")



    x_vals = []
    y_vals = []
    for measurement in All_measurements["ArconF"]:
        if 6 < measurement.pH < 7:
            x_vals.append(measurement.TgT)
            y_vals.append(measurement.Plateau)
            plt.errorbar(measurement.TgT, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"bo")

    # params = fit_equation(x_vals, y_vals, True)
    # points = linear_eq(x_indices, *params)
    # plt.plot(x_indices, 10**points, "b--", label="fit", zorder=10)

    ax.set_yscale("log")

    plt.ylabel("G0 [kPa]")
    plt.xlabel("Tg/T [-]")
    ax.legend(handles=[Patch(facecolor='black',
                             label='Alpha 8'),
                       Patch(facecolor='blue', label='Arcon F')])
    plt.title("TgT")
    plt.ylim([1e0, 1e3])
    plt.show()





def pH_comparison():

    fig, ax1 = plt.subplots(1, 1, dpi=600)
    # for measurement in measurements:
    #    if measurement.Temperature==170.0 and measurement.Moisture == 0.55:
    #        plt.errorbar(measurement.pH, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"bo")
    for measurement in All_measurements["ArconF"]:
        if measurement.Temperature == 150.0 and measurement.Moisture == 0.55 and not 6.3 < measurement.pH <6.6:
            ax1.errorbar(measurement.pH, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"bo")

    ax1.set_yscale("log")
    ax1.set_ylabel("G0 [kPa]")
    ax1.set_xlabel("pH")


    for i, measurement in enumerate(All_measurements["Alpha8"]):
        if measurement and measurement.Temperature == 150.0 and measurement.Moisture == 0.55:
            print(i, measurement[:-1])
            ax1.errorbar(measurement.pH, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"ko")


    ax1.set_title("pH comparison")
    ax1.legend(handles=[Patch(facecolor='black',label='Alpha 8'),Patch(facecolor='blue', label='Arcon F')])
    plt.ylim([4e0, 3e1])
    plt.show()


def Plateau_comparison():
    x_indices = np.linspace(90, 170, 20)

    x_vals = []
    y_vals = []

    fig, ax = plt.subplots()
    for measurement in All_measurements["ArconF"]:
        if 6 < measurement.pH < 7 and measurement.Moisture == 0.55:
            x_vals.append(measurement.Temperature)
            y_vals.append(measurement.Plateau)
            plt.errorbar(measurement.Temperature, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"bo")
    #params = fit_equation(x_vals, y_vals, True)
    #points = linear_eq(x_indices, *params)
    #plt.plot(x_indices, 10**points, "b--", label="fit", zorder=10)


    x_vals = []
    y_vals = []
    for measurement in All_measurements["Alpha8"]:
        if measurement and 6 < measurement.pH < 7  and measurement.Moisture == 0.55:
            x_vals.append(measurement.Temperature)
            y_vals.append(measurement.Plateau)
            plt.errorbar(measurement.Temperature, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"ko")
    #params = fit_equation(x_vals, y_vals, True)
    #points = linear_eq(x_indices, *params)
    #plt.plot(x_indices, 10**points, "k--", label="fit", zorder=10)
    ax.legend(handles=[Patch(facecolor='black',
                             label='Alpha 8'),
                       Patch(facecolor='blue', label='Arcon F')])
    ax.set_yscale("log")


    plt.ylabel("G0 [kPa] ")
    plt.xlabel("T [°C]")
    plt.title("plateau")
    plt.show()

def tan_delta_master():
    # fig, ax = plt.subplots()
    # for measurement in All_measurements["Alpha8"]:
    #     if measurement and 6 < measurement.pH < 7  and measurement.Moisture == 0.55:
    #         plt.errorbar(measurement.Temperature, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"ko")
    #
    # ax.set_yscale("log")
    # plt.ylabel("G0 [kPa] ")
    # plt.xlabel("T [°C]")
    # plt.title("plateau")
    # plt.show()
    #
    # plt.figure()
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, layout="constrained")
    for measurement in All_measurements["Alpha8"]:
        if measurement and 6 < measurement.pH < 7  and measurement.Moisture == 0.55:
            ax1.loglog(measurement.data["Strain %"], measurement.data["G* [kPa]"], c=plt.cm.viridis((measurement.Temperature-80)/90), label=measurement.Info)
    # for measurement in All_measurements["Alpha8"]:
    #     if measurement and 6 < measurement.pH < 7  and measurement.Moisture == 0.55:
    #         ax1.loglog(measurement.data["Strain %"], measurement.data["G' [kPa]"]/measurement.Plateau,
    #                    c=plt.cm.viridis((measurement.Temperature-80)/90), label=measurement.Info)
    #         ax1.loglog(measurement.data["Strain %"], measurement.data["G'' [kPa]"] / measurement.Plateau,
    #                    c=plt.cm.viridis((measurement.Temperature - 80) / 90), label=measurement.Info)
    #
    for measurement in All_measurements["ArconF"]:
        if measurement and 6 < measurement.pH < 7  and measurement.Moisture == 0.55:
                ax2.loglog(measurement.data["Strain %"], measurement.data["G* [kPa]"], c=plt.cm.viridis((measurement.Temperature-80)/90), label=measurement.Info)
    # for measurement in All_measurements["ArconF"]:
    #     if measurement and 6 < measurement.pH < 7  and measurement.Moisture == 0.55:
    #         ax2.loglog(measurement.data["Strain %"], measurement.data["G' [kPa]"]/measurement.Plateau,
    #                    c=plt.cm.viridis((measurement.Temperature-80)/90), label=measurement.Info)
    #         ax2.loglog(measurement.data["Strain %"], measurement.data["G'' [kPa]"] / measurement.Plateau,
    #                    c=plt.cm.viridis((measurement.Temperature - 80) / 90), label=measurement.Info)
    fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=80, vmax=170), cmap=plt.cm.viridis),
                 ax=ax2, orientation='vertical', label='T (°)')

    ax1.set_title("Alpha 8")
    ax2.set_title("Arcon F")
    #ax1.set_ylabel("G'/G0 G\"/G0[kPa]")
    ax1.set_ylabel("G* [kPa]")
    ax1.set_xlabel("strain [%]")
    ax2.set_xlabel("strain [%]")
    ax1.set_ylim(1e-1, 5e2)
    plt.show()

print(Measurement._fields)
Drive = "D"
ArconF_T()
#new_Alpha8()
#pH_comparison()
#TgT_comparison()
#Plateau_comparison()
JMP_export()
#tan_delta_master()

