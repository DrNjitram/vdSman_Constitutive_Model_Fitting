import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import namedtuple

from matplotlib.patches import Patch
from scipy.optimize import curve_fit

Measurement = namedtuple("Measurement", "Mixture, Customer, Instrument, Material, Comments, Info, Batch_NB, Order_NB, N3, N4, N5, Frequency, Temperature, InterID, Moisture, pH, TgT, Plateau, Plateau_std, data")
#path = r"D:\Data\Rheology\CCR\202409 DOE pH T M Arcon F\Overview.xlsx"
colors = ["k", "r", "g", "b", "m", "c", "y", "k", "k", "k", "k", "k", "k"]
Cpw = 1.92E+03
Tgw = 139 #K
Cps = 4.25E+02
Tgs = 384 #K
debug = False
All_measurements = dict()
unique_id = []
fit_eq = lambda x, *params: params[0] + params[1]*(1-np.exp(-x))

# pH correction data
# 2nd order poly
base_data = [1.7891, 3.9621, 6.3907]
acid_data = [4.2209, -5754, 6.4553]

def fit_equation(x_vals, y_vals):
    params = [0, 1]
    params, covariance = curve_fit(fit_eq, x_vals, np.log10(y_vals), p0=params)
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


def get_measurements(path, filenames, measurement_width, header_row, headers, default_pH=7.0, re_calculate_pH=False):
    measurements = []
    for filename in filenames:
        index = 0
        file_path = os.path.join(path, filename)
        print(filename)
        while True:
            try:
                info = pd.read_excel(io=file_path, header=None, skiprows=1, nrows=14, usecols=[(1 if not index else 0)+index*measurement_width], names=["Value"])
                data = pd.read_excel(io=file_path, header=None, names=headers, skiprows=header_row, usecols=list(range(index * measurement_width, (index+1)*measurement_width)))
                measurement = Measurement(*info["Value"], -1, default_pH, -1, -1, -1, data)
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
                                pH = bit
                            elif "NaOH" in information_string:
                                pH = -1*bit
                            elif bit != 0.0:
                                if bit == 0.752:
                                    pH = -0.52
                                elif bit in [0.07, 0.24, 0.41, 0.74]:
                                    pH = bit
                                else:
                                    pH = -1*bit
                            else:
                                pH = 0.0

                            measurement = measurement._replace(pH=float(pH))


                    measurement = measurement._replace(Temperature=float(measurement.Temperature[:-2]))
                    measurement = measurement._replace(Plateau=np.mean(measurement.data["G' [kPa]"][5:10]))
                    measurement = measurement._replace(Plateau_std=np.std(measurement.data["G' [kPa]"][5:10]))


                    if "T Run" in filename: # Alpha 8
                        measurement = measurement._replace(TgT=float(measurement.Material))
                        measurement = measurement._replace(Moisture=backcalculate_moisture(measurement.TgT, measurement.Temperature))
                    elif r"pH.xlsx" in filename:
                        measurement = measurement._replace(pH=recalculate_pH(measurement.Material))
                        measurement = measurement._replace(Moisture=0.6)
                        measurement = measurement._replace(TgT=TgT(measurement.Moisture, measurement.Temperature))
                    else:
                        measurement = measurement._replace(TgT = TgT(measurement.Moisture, measurement.Temperature))

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
    measurement_width = 19
    header_row = 18
    headers = ['Strain %', "S' [dNm]", "S'' [dNm]", 'S* [dNm]', "G' [kPa]", "G'' [kPa]", 'G* [kPa]', 'Tan Delta',
               'Temperature [°C]', 'Pressure [bar]', "Viscosity n' [Pa*s]", 'Viscosity n" [Pa*s]',
               'Viscosity n* [Pa*s]', 'Shear Rate [1/sec]', 'I2 abs.', 'I3 abs.', 'I2/I1 [%]', 'I3/I1 [%]',
               'Shear Stress [Pa]']

    path = rf"{Drive}:\Data\Rheology\CCR\202409 DOE pH T M Arcon F"
    filenames = ["Overview_fixed.xls", "Overview2.xls", "Overview3.xls", "Overview4.xls", ]
    measurements = get_measurements(path, filenames, measurement_width, header_row, headers, 6.42, True)
    measurements.pop(11)


    measurement_width = 12
    headers = ['Strain %', "G' [kPa]", "G'' [kPa]", 'G* [kPa]', 'Tan Delta',
               'Temperature [°C]', 'Pressure [bar]', "Viscosity n' [Pa*s]", 'Viscosity n" [Pa*s]',
               'Viscosity n* [Pa*s]', 'Shear Rate [1/sec]', 'Shear Stress [Pa]']
    measurements2 = get_measurements(path, ["Overview5.xlsx", "Overview6.xlsx"], measurement_width, header_row, headers, 6.42, False)
    measured_pH = [4.80, 5.01, 5.87, 6.02, 4.78, 5.61, 5.22, 5.17, 5.60, 6.31, 6.46, 6.45, 6.42, 6.85 ,8.62, 9, 5.03, 6.26]
    for measurement, pH in zip(measurements2[1:], measured_pH[1:]): # First one is busted
        #print(f"{pH}, {measurement.Comments}, {measurement.Info}, {measurement.Temperature}")
        measurement = measurement._replace(pH=pH)
        measurements.append(measurement)

    All_measurements["ArconF"] = measurements

    if debug or True: # To export into JMP
        print("T (K), M (%), pH, TgT (-), G0 [kPa], log10 G0")
        for measurement in measurements:
            print(f"{measurement.Temperature+273.15}, {measurement.Moisture:.2f}, {measurement.pH:.2f}, {measurement.TgT:.2f}, {measurement.Plateau:.2f}, {np.log10(measurement.Plateau):.2f}")

    if debug:
        plt.figure()
        X = "Strain %"
        i = 0
        for measurement in measurements:
            if measurement.Moisture == 0.55 and 6.5 < measurement.pH < 7.5:
                plt.loglog(measurement.data[X], measurement.data["G' [kPa]"], f"{colors[i]}o", label=measurement.Temperature)
                plt.loglog(measurement.data[X], measurement.data["G'' [kPa]"], f"{colors[i]}s")
                i += 1

        plt.ylim([10**-2, 10**2])
        plt.ylabel("G' [kPa] / G'' [kPa]")
        plt.xlabel(X)
        plt.legend()
        plt.title("M55% pH7")
        plt.show()

    plt.figure()

    graph_measurement(measurements[34])

    graph_measurement(measurements[27])

    fig, ax = plt.subplots()
    for measurement in measurements:
        if  6 < measurement.pH < 7 and measurement.Moisture == 0.55:
            plt.errorbar(measurement.Temperature, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"bo")
    ax.set_yscale("log")
    plt.ylabel("G0 [kPa] ")
    plt.xlabel("T [°C]")
    plt.title("plateau")
    plt.show()

    fig, ax = plt.subplots()
    #for measurement in measurements:
    #    if measurement.Temperature==170.0 and measurement.Moisture == 0.55:
    #        plt.errorbar(measurement.pH, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"bo")
    for measurement in measurements:
        if measurement.Temperature==150.0 and measurement.Moisture == 0.55:
            plt.errorbar(measurement.pH, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"ro")

    ax.set_yscale("log")
    plt.ylabel("G0 [kPa] ")
    plt.xlabel("pH")
    plt.title("Arcon F pH comparison")
    #ax.legend(handles=[Patch(facecolor='blue',label='T=170'),Patch(facecolor='red', label='T=150')])
    plt.show()

    fig, ax = plt.subplots()
    # for measurement in measurements:
    #    if measurement.Temperature==170.0 and measurement.Moisture == 0.55:
    #        plt.errorbar(measurement.pH, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"bo")
    i = 0
    for idx in [37, 39, 40, 45]:
        measurement = measurements[idx]
        if measurement.Temperature == 150.0 and measurement.Moisture == 0.55 and measurement.pH < 6.5:
            norm = measurement.Plateau
            plt.loglog(measurement.data["Strain %"], measurement.data["G' [kPa]"]/norm, f"{colors[i]}o", label=measurement.pH)
            plt.loglog(measurement.data["Strain %"], measurement.data["G'' [kPa]"]/norm, f"{colors[i]}s")
            i += 1

    plt.ylim([10**-3, 1])
    plt.xlim([10**0, 10**2.5])
    plt.ylabel("normalised G' [kPa] ")
    plt.xlabel("Strain %")
    plt.title("Arcon F pH comparison 2")
    plt.legend()
    # ax.legend(handles=[Patch(facecolor='blue',label='T=170'),Patch(facecolor='red', label='T=150')])
    plt.show()

def Alpha8_T():
    measurement_width = 22
    header_row = 18
    headers = ['Strain %', "S' [dNm]", "S'' [dNm]", 'S* [dNm]', "G' [kPa]", "G'' [kPa]", 'G* [kPa]', 'Tan Delta', 'Temperature [°C]', 'Pressure [bar]', "Viscosity n' [Pa*s]", 'Viscosity n" [Pa*s]', 'Viscosity n* [Pa*s]', 'Shear Rate [1/sec]', 'J\' [1/kPa]', 'J" [1/kPa]','J* [1/kPa]', 'I2 abs.', 'I3 abs.', 'I2/I1 [%]', 'I3/I1 [%]', 'Shear Stress [Pa]']
    path = rf"{Drive}:\Data\Rheology\CCR\2023-05-19 pH and T run\T"
    filenames = ["T Run.xls"]
    measurements = get_measurements(path, filenames, measurement_width, header_row, headers, 7.0)
    All_measurements["Alpha8_T"] = measurements

    measurement = measurements[0]
    plt.loglog(measurement.data["Strain %"], measurement.data["G' [kPa]"], f"{colors[0]}o", label="G'")
    plt.loglog(measurement.data["Strain %"], measurement.data["G'' [kPa]"], f"{colors[0]}s", label="G\"")
    plt.ylabel("G' [kPa]/G\" [kPa] ")
    plt.xlabel("Strain %")
    plt.legend()
    plt.title(
        f"T{measurement.Temperature}°C, M{measurement.Moisture * 100:.0f}%, pH{measurement.pH:.2f}")
    plt.show()

    plt.figure()
    X = "Strain %"
    for i, idx in enumerate(range(len(measurements))):
        measurement = measurements[idx - 1]
        plt.loglog(measurement.data[X], measurement.data["G' [kPa]"], f"{colors[i]}o", label=measurement.Material)
        plt.loglog(measurement.data[X], measurement.data["G'' [kPa]"], f"{colors[i]}s")


    plt.ylabel("G' [kPa] / G'' [kPa]")
    plt.xlabel(X)
    plt.title("Alpha 8 T")
    plt.legend()
    plt.show()

def old_Alpha8_pH():
    measurement_width = 22
    header_row = 18
    headers = ['Strain %', "S' [dNm]", "S'' [dNm]", 'S* [dNm]', "G' [kPa]", "G'' [kPa]", 'G* [kPa]', 'Tan Delta', 'Temperature [°C]', 'Pressure [bar]', "Viscosity n' [Pa*s]", 'Viscosity n" [Pa*s]', 'Viscosity n* [Pa*s]', 'Shear Rate [1/sec]', 'J\' [1/kPa]', 'J" [1/kPa]','J* [1/kPa]', 'I2 abs.', 'I3 abs.', 'I2/I1 [%]', 'I3/I1 [%]', 'Shear Stress [Pa]']
    path = rf"{Drive}:\Data\Rheology\CCR\2023-05-19 pH and T run\pH"
    filenames = ["pH.xlsx"]
    measurements = get_measurements(path, filenames, measurement_width, header_row, headers, True)

    All_measurements["Alpha8_pH"] = measurements

    for measurement in measurements:
        plt.semilogy(float(measurement.pH), measurement.Plateau, f"ko")

    plt.title("Alpha 8 pH")
    plt.ylabel("plateau G' [kPa] ")
    plt.xlabel("pH")
    plt.show()

def TgT_comparison():
    fig, ax = plt.subplots()
    x_indices = np.linspace(0.35, 0.55, 20)

    x_vals = []
    y_vals = []
    for measurement in All_measurements["Alpha8"]:
        if measurement and measurement.pH == 0.0:
            x_vals.append(measurement.TgT)
            y_vals.append(measurement.Plateau)
            plt.errorbar(float(measurement.TgT), measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"ko")

    params = fit_equation(x_vals, y_vals)
    points = fit_eq(x_indices, *params)
    plt.plot(x_indices, 10**points, "k--", label="fit", zorder=10)

    x_vals = []
    y_vals = []
    for measurement in All_measurements["ArconF"]:
        if 6 < measurement.pH < 7:
            x_vals.append(measurement.TgT)
            y_vals.append(measurement.Plateau)
            plt.errorbar(float(measurement.TgT), measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"bo")

    params = fit_equation(x_vals, y_vals)
    points = fit_eq(x_indices, *params)
    plt.plot(x_indices, 10**points, "b--", label="fit", zorder=10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.ylabel("plateau G' [kPa] ")
    plt.xlabel("Tg/T [-]")
    ax.legend(handles=[Patch(facecolor='black',
                             label='Alpha 8'),
                       Patch(facecolor='blue', label='Arcon F')])
    plt.title("TgT")
    plt.show()

def new_Alpha8():
    measurement_width = 12
    header_row = 18
    headers = ['Strain %', "G' [kPa]", "G'' [kPa]", 'G* [kPa]', 'Tan Delta', 'Temperature [°C]', 'Pressure [bar]', "Viscosity n' [Pa*s]", 'Viscosity n" [Pa*s]', 'Viscosity n* [Pa*s]', 'Shear Rate [1/sec]', 'Shear Stress [Pa]']
    path = rf"{Drive}:\Data\Rheology\CCR\202411 DoE pH T M Alph 8"
    filenames = ["Overview 1.xls", "Overview 2.xls", "Overview 3.xls", "Overview 4.xls"]
    measurements = get_measurements(path, filenames, measurement_width, header_row, headers, default_pH=0)
    measurements[7] = False
    measurements[17] = False
    measurements[25] = False
    All_measurements["Alpha8"] = measurements

    fig, ax = plt.subplots()
    for measurement in measurements:
        if  measurement and measurement.pH == 0.0:
            plt.errorbar(float(measurement.TgT), measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"ko")

    ax.set_yscale("log")
    plt.title("Alpha 8")
    plt.ylabel("plateau G' [kPa] ")
    plt.xlabel("Tg/T")
    plt.show()

    fig, ax = plt.subplots()
    for i, measurement in enumerate(measurements):
        if measurement and measurement.Temperature == 150.0 and measurement.Moisture ==0.55:
            print(i, measurement[:-1])
            plt.errorbar(measurement.pH, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"ko")

    ax.set_yscale("log")
    plt.title("Alpha 8")
    plt.ylabel("plateau G' [kPa] ")
    plt.xlabel("pH")
    plt.show()

def pH_comparison():
    measurements = All_measurements["ArconF"]
    fig, ax = plt.subplots()
    # for measurement in measurements:
    #    if measurement.Temperature==170.0 and measurement.Moisture == 0.55:
    #        plt.errorbar(measurement.pH, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"bo")
    for measurement in measurements:
        if measurement.Temperature == 150.0 and measurement.Moisture == 0.55:
            plt.errorbar(measurement.pH, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"ro")

    ax.set_yscale("log")
    plt.ylabel("G0 [kPa] ")
    plt.xlabel("pH")
    plt.title("Arcon F pH")

    measurements = All_measurements["Alpha8"]
    fig, ax = plt.subplots()
    for i, measurement in enumerate(measurements):
        if measurement and measurement.Temperature == 150.0 and measurement.Moisture == 0.55:
            print(i, measurement[:-1])
            plt.errorbar(measurement.pH, measurement.Plateau, yerr=measurement.Plateau_std, fmt=f"ko")

    ax.set_yscale("log")
    plt.title("Alpha 8")
    plt.ylabel("plateau G' [kPa] ")
    plt.xlabel("Alpha 8 pH")
    plt.show()





print(Measurement._fields)
Drive = "W"
ArconF_T()
#Alpha8_T()
#old_Alpha8_pH()

new_Alpha8()
pH_comparison()
TgT_comparison()

