# (C) 2025 Francesco Paparella & Faiq Raees
# This file is part of Sea_Ice_RAP.
#
# Sea_Ice_RAP is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Sea_Ice_RAP is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Sea_Ice_RAP. If not, see <https://www.gnu.org/licenses/>.

from numpy import *
from scipy.special import factorial
import warnings

def lowess(x, y, x_resample=None,
           deltax=None, degree=2, weight_function=None, niter=2):
    """Robust locally weighted regression and smoothing scatterplots.
See: W. S. Cleveland Journal of the American Statistical Association,
Vol. 74, No. 368. (1979), 829-836.

The main difference with the algorithm described in the paper is that
each local regression occurs on a window of fixed width 'deltax',
rather than on a fixed fraction 'f' of the data points. This is
preferable when the data is unevenly sampled, and one desires to
smooth out the data on scales much smaller than 'deltax', while
retaining a good resolution on the fluctuations at that or larger
scales. Because unevenly sampled data are assumed, this implementation
performs a complete least-square fit for each local regression, and is
therefore slower than implementations that exploit regular sampling.

    """
    def check_enough_points_in_window(idx):
        if len(idx) < degree+1:
            raise ValueError("The window deltax is too narrow: not enough data to fit the required polynomial!\nYou need at least (degree+1) points in your window, and at the beginning and the end of the time series the window length is just deltax/2. Recall that window edges are excluded (their weight is zero, anyway).\n(Also check for holes in the data)")

    def pade(x):
        return 1./(1+2*(x*x)+3*(x*x*x*x))

    def bisquare(x):
        a = (1. - x*x)
        a *= a>0
        return a*a

    def tricube(x):
        a = (1. - x*x*abs(x))
        a *= a>0
        return a*a*a

    if deltax==None:
        deltax = (amax(x)-amin(x))/4.
    if weight_function==None:
        weight_function = tricube
    delta = ones_like(x)
    dhalf = deltax/2.
    for nloop in range(niter):
        e = zeros_like(x)
        for i, xi in enumerate(x):
            idx = where((x > xi-dhalf)*(x < xi+dhalf))[0]
            check_enough_points_in_window(idx)
            xfit = x[idx] - xi
            yfit = y[idx]
            wfit = delta[idx] * weight_function(xfit/dhalf)

            poly = polyfit(xfit, yfit, degree, w=wfit)
            e[i] = y[i] - poly[-1]
        delta = pade(e/(6*median(abs(e))))


    if x_resample is None:
        x_resample = x
    derivatives = []
    cfc = factorial(arange(degree+1))
    for xi in x_resample:
        idx = where((x > xi-dhalf)*(x < xi+dhalf))[0]
        check_enough_points_in_window(idx)
        xfit = x[idx] - xi
        yfit = y[idx]
        wfit = delta[idx] * weight_function(xfit/dhalf)
        poly = polyfit(xfit, yfit, degree, w=wfit)
        derivatives.append(poly[::-1]*cfc)
    return array(derivatives)

