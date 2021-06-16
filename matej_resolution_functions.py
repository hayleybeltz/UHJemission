import numpy as np
import sys
import astropy.constants as ac


K_B = ac.k_B.cgs.value
H = ac.h.cgs.value
C = ac.c.cgs.value


def percent_counter(z, nz, y=0, ny=1, x=0, nx=1):
    """ displays percentage completed of a long operation (usually a for loop) for up to three indices """

    percentage = float((x + nx * y + nx * ny * z) / (nx * ny * nz) * 100.0)
    sys.stdout.write("calculating: {:.1f}%\r".format(percentage))
    sys.stdout.flush()


def calc_analyt_planck_in_interval(temp, lower_lambda, higher_lambda):
    """ calculates the planck function over a wavelength integral

    :param temp: float
                 the blackbody temperature

    :param lower_lambda: float
                         lower wavelength boundary (cm units!)

    :param higher_lambda: float
                          upper wavelength boundary (cm units!)

    :return: float
             the Planckian blackbody function integrated over a wavelength interval and averaged
    """

    d = 2.0 * (K_B / H)**3 * K_B * temp**4 / C**2
    y_top = H * C / (higher_lambda * K_B * temp)
    y_bot = H * C / (lower_lambda * K_B * temp)

    result = 0

    for n in range(1, 200):  # 200 found to be accurate enough.
        result += np.exp(-n*y_top) * (y_top**3/n + 3.0*y_top**2/n**2 + 6.0*y_top/n**3 + 6.0/n**4) \
                - np.exp(-n*y_bot) * (y_bot**3/n + 3.0*y_bot**2/n**2 + 6.0*y_bot/n**3 + 6.0/n**4)

    result *= d / (higher_lambda - lower_lambda)

    return result


def gauss_pdf(x, x_0, hwhm):

    pdf = (np.log(2)/np.pi)**0.5 / hwhm * np.exp(-np.log(2) * ((x - x_0) / hwhm)**2)

    return pdf


def convolve_with_gaussian(old_lamda, old_flux, resolution, new_lamda=None):

    if new_lamda is None:

        new_lamda = [old_lamda[0]]

        while new_lamda[-1] < old_lamda[-1]:

            new_lamda.append(new_lamda[-1] + new_lamda[-1] / resolution)

    # need wavelength widths of old lambda grid
    delta_lamda = np.zeros(len(old_lamda))

    for ll in range(len(old_lamda)):

        if ll == 0:

            delta_lamda[ll] = (old_lamda[ll + 1] - old_lamda[ll])

        elif ll == len(old_lamda) - 1:

            delta_lamda[ll] = (old_lamda[ll] - old_lamda[ll - 1])

        else:
            delta_lamda[ll] = (old_lamda[ll + 1] - old_lamda[ll - 1]) / 2

    # calculation of new flux values due to convolution with older ones
    flux_conv = np.zeros(len(new_lamda))

    for l in range(len(new_lamda)):

        percent_counter(l, len(new_lamda))

        # FWHM of Gaussian pdf equals the resolving power R (and thus HWHM = R/2)
        hwhm = new_lamda[l] / (2 * resolution)

        for ll in range(len(old_lamda)):

            if old_lamda[ll] - new_lamda[l] < -5 * hwhm:
                continue

            elif old_lamda[ll] - new_lamda[l] > 5 * hwhm:
                break

            else:
                flux_conv[l] += old_flux[ll] * gauss_pdf(new_lamda[l], old_lamda[ll], hwhm) * delta_lamda[ll]

    return new_lamda, flux_conv


def convert_spectrum(old_lambda, old_flux, new_lambda, int_lambda=None, type='log', extrapolate_with_BB_T=0):
    """ converts a spectrum from one to another resolution. This method conserves energy. It is the real deal.

        :param old_lambda: list of float or numpy array
                           wavelength values of the old grid to be discarded. must be in ascending order!

        :param old_flux: list of float or numpy array
                         flux values of the old grid to be discarded.

        :param new_lambda: list of float or numpy array
                           wavelength values of the new output grid. must be in ascending order!

        :param int_lambda: (optional) list of float or numpy array
                           wavelength values of the interfaces of the new grid bins. must be in ascending order!
                           if not provided they are calculated by taking the middle points between the new_lambda values.

        :param type: (optional) 'linear' or 'log'
                     either linear interpolation or logarithmtic interpolation possible

        :param extrapolate_with_BB_T: (optional) float
                                      the out-of-boundary flux values will be extrapolated with a blackbody spectrum.
                                      set here the temperature. if not provided the out-of-boundary flux values are set to zero.

        :return: list of floats
                 flux values at the new wavelength grid points

        """

    if int_lambda is None:

        int_lambda = []

        int_lambda.append(new_lambda[0] - (new_lambda[1] - new_lambda[0]) / 2)

        for x in range(len(new_lambda) - 1):
            int_lambda.append((new_lambda[x + 1] + new_lambda[x]) / 2)

        int_lambda.append(new_lambda[-1] + (new_lambda[-1] - new_lambda[-2]) / 2)

    if extrapolate_with_BB_T > 0:

        extrapol_values = []

        #print("\n\nPre-tabulating blackbody values with a temperature of " + str(extrapolate_with_BB_T) + " K ...\n")

        for i in range(len(new_lambda)):

            percent_counter(i, len(new_lambda))

            extrapol_values.append(np.pi * calc_analyt_planck_in_interval(extrapolate_with_BB_T, int_lambda[i], int_lambda[i+1]))

        #print("Pre-tabulation done! \n")

    elif extrapolate_with_BB_T == 0:

        extrapol_values = np.zeros(len(new_lambda))

    else:
        raise ValueError("Error: extrapolation blackbody temperature cannot be negative.")

    #print("Starting pre-conversion...\n")

    int_flux = [0] * len(int_lambda)
    new_flux = []
    old_lambda = np.array(old_lambda)  # conversion required by np.where

    if type == 'linear':

        for i in range(len(int_lambda)):

            percent_counter(i, len(int_lambda))

            if int_lambda[i] < old_lambda[0]:
                continue

            elif int_lambda[i] > old_lambda[len(old_lambda) - 1]:
                break

            else:
                p_bot = len(np.where(old_lambda < int_lambda[i])[0]) - 1

                interpol = old_flux[p_bot] * (old_lambda[p_bot + 1] - int_lambda[i]) + old_flux[p_bot + 1] * (int_lambda[i] - old_lambda[p_bot])
                interpol /= (old_lambda[p_bot + 1] - old_lambda[p_bot])
                int_flux[i] = interpol

        #print("\n  Pre-conversion done!")

        #print("\nStarting main conversion...\n")

        for i in range(len(new_lambda)):

            percent_counter(i, len(new_lambda))

            if int_flux[i] == 0 or int_flux[i+1] == 0:
                new_flux.append(extrapol_values[i])

            else:
                p_bot = len(np.where(old_lambda < int_lambda[i])[0]) - 1

                p_start = p_bot + 1

                for p in range(p_start, len(old_lambda)):

                    if p == p_start:
                        if old_lambda[p_start] < int_lambda[i + 1]:
                            interpol = (int_flux[i] + old_flux[p]) / 2.0 * (old_lambda[p] - int_lambda[i])

                        else:
                            interpol = (int_flux[i] + int_flux[i + 1]) / 2.0
                            break
                    else:
                        if old_lambda[p] < int_lambda[i + 1]:
                            interpol += (old_flux[p - 1] + old_flux[p]) / 2.0 * (old_lambda[p] - old_lambda[p - 1])

                        else:
                            interpol += (old_flux[p - 1] + int_flux[i + 1]) / 2.0 * (int_lambda[i + 1] - old_lambda[p - 1])
                            interpol /= (int_lambda[i + 1] - int_lambda[i])
                            break

                new_flux.append(interpol)

    elif type == 'log':

        for i in range(len(int_lambda)):

            percent_counter(i, len(int_lambda))

            if int_lambda[i] < old_lambda[0]:
                continue

            elif int_lambda[i] > old_lambda[len(old_lambda) - 1]:
                break

            else:
                p_bot = len(np.where(old_lambda < int_lambda[i])[0]) - 1

                interpol = old_flux[p_bot] ** (old_lambda[p_bot + 1] - int_lambda[i]) * old_flux[p_bot + 1] ** (int_lambda[i] - old_lambda[p_bot])
                interpol = interpol**(1/(old_lambda[p_bot + 1] - old_lambda[p_bot]))
                int_flux[i] = interpol

        #print("\n  Pre-conversion done!")

        #print("\nStarting main conversion...\n")

        for i in range(len(new_lambda)):

            percent_counter(i, len(new_lambda))

            if int_flux[i] == 0 or int_flux[i + 1] == 0:
                new_flux.append(extrapol_values[i])

            else:
                p_bot = len(np.where(old_lambda < int_lambda[i])[0]) - 1

                p_start = p_bot + 1

                for p in range(p_start, len(old_lambda)):

                    if p == p_start:
                        if old_lambda[p_start] < int_lambda[i + 1]:
                            interpol = (int_flux[i] * old_flux[p]) ** (0.5 * (old_lambda[p] - int_lambda[i]))

                        else:
                            interpol = (int_flux[i] * int_flux[i + 1]) ** 0.5
                            break
                    else:
                        if old_lambda[p] < int_lambda[i + 1]:
                            interpol *= (old_flux[p - 1] * old_flux[p]) ** (0.5 * (old_lambda[p] - old_lambda[p - 1]))

                        else:
                            interpol *= (old_flux[p - 1] * int_flux[i + 1]) ** (0.5 * (int_lambda[i + 1] - old_lambda[p - 1]))
                            interpol = interpol ** (1/(int_lambda[i + 1] - int_lambda[i]))
                            break

                new_flux.append(interpol)

    #print("\n  Main conversion done!")

    return new_flux


def rebin_spectrum_to_resolution(old_lamda, old_flux, resolution, w_unit='cm', type='log'):
    """ rebins a given spectrum to a new resolution

    :param old_lambda:  list of float or numpy array
                           wavelength values of the old grid to be discarded.
                           must be in ascending order!

    :param old_flux:    list of float or numpy array
                        flux values of the old grid to be discarded.

    :param resolution:  float
                        resolution (R=lamda/delta_lamda) of the new, rebinned wavelength grid

    :param w_unit:      'cm'/'micron'
                        the units of the old wavelength values. will be the output units as well.

    :param type:        (optional) 'linear', 'log' or 'gaussian'
                        - linear interpolation is the default. it conserves the total energy in each bin.
                        - logarithmic interpolation is usually used for opacities or other quantities where the integral does not need to be conserved.
                        - alternatively, the spectrum can be convolved with a Gaussian distribution, where FWHM = R

    :return 1:          wavelength values of the rebinned grid
    :return 2:          flux values of the rebinned grid
    """

    if w_unit == 'micron':
        old_lamda = [l * 1e-4 for l in old_lamda]

    bot_limit = old_lamda[0]

    top_limit = old_lamda[-1]

    rebin_lamda = []
    l_point = bot_limit

    while l_point < top_limit:
        rebin_lamda.append(l_point)

        l_point *= (resolution + 1) / resolution

    if type == "gaussian":
        _, rebin_flux = convolve_with_gaussian(old_lamda, old_flux, resolution, rebin_lamda)
    else:
        rebin_flux = convert_spectrum(old_lamda, old_flux, rebin_lamda, type=type, extrapolate_with_BB_T=0)

    if w_unit == 'micron':
        rebin_lamda = [l * 1e4 for l in rebin_lamda]

    return rebin_lamda, rebin_flux
