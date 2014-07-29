import numpy as np
import obspy.core.trace as trace
from scipy.fftpack import fft, ifft

def calcWindow(waveform, windowStart, opt):

    """
    Calculates the amplitude coefficient and FFT for a window of data.

    waveform: numpy array of waveform data
    windowStart: starting sample of window
    opt: Options object describing station/run parameters

    Returns windowCoeff and windowFFT
    """
    
    # Calculate filtered version of data for use in cross correlation
    trig = trace.Trace(waveform, {"sampling_rate":opt.samprate})
    waveform = trig.filter("bandpass", freqmin=opt.fmin, freqmax=opt.fmax, corners=2,
		zerophase=True).data

    windowCoeff = 1/np.sqrt(sum(waveform[windowStart:(windowStart + opt.winlen)] *
        waveform[windowStart:(windowStart + opt.winlen)]))
    windowFFT = np.reshape(fft(waveform[windowStart:(windowStart + opt.winlen)]),
        (opt.winlen,))

    return windowCoeff, windowFFT


def xcorr1x1(windowFFT1, windowFFT2, windowCoeff1, windowCoeff2):

    """
    Calculates the cross-correlation coefficient and lag for two windows.

    windowFFT1: FFT of first window
    windowFFT2: FFT of second window
    windowCoeff1: amplitude coefficient of first window
    windowCoeff2: amplitude coefficient of second window

    Order matters for sign of lag, but not CCC.

    Returns maximum cross-correlation and optimal lag (in samples)
    """

    M = len(windowFFT1)
    coeff = windowCoeff1 * windowCoeff2

    lags = np.roll(np.linspace(-M/2 + 1, M/2, M, endpoint=True),
        M/2 + 1).astype(int)
    cors = np.real(ifft(windowFFT1 * np.conj(windowFFT2))) * coeff

    indx = np.argmax(cors)
    maxcor = cors[indx]
    maxlag = lags[indx]

    return maxcor, maxlag