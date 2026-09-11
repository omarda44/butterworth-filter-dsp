# =============================================================
#  Project #13 - Butterworth filter of order N, applied to x[n]
#                entirely in the frequency domain.
#
#  Method:
#    1. Prewarp the cutoff, place the analog Butterworth poles.
#    2. Map s -> z with the bilinear transform, evaluated on the
#       unit circle to get H[k] at the DFT bin frequencies.
#    3. Y[k] = H[k] * X[k], then y[n] = IDFT{Y[k]}.
# =============================================================

import numpy as np
import matplotlib.pyplot as plt

SAMPLING_FREQUENCY = 1000.0     # fs [Hz]
CUTOFF_FREQUENCY = 100.0        # fc [Hz]


def ask(message, parser):
    """Keep asking until the user types something valid."""
    while True:
        try:
            return parser(input(message))
        except (ValueError, SyntaxError, IndexError):
            print("  invalid input, try again.")


def parse_order(text):
    order = int(text)
    if not 1 <= order <= 20:
        raise ValueError
    return order


def parse_signal(text):
    values = [float(v) for v in text.strip(" []()").split(",") if v.strip()]
    if len(values) < 1:
        raise ValueError
    return np.array(values)


# ---------------- 1. Input ----------------
filter_order = ask("Enter N: ", parse_order)
input_signal = ask("Enter x[n]: ", parse_signal)

number_of_samples = len(input_signal)
sample_indices = np.arange(number_of_samples)

# The DFT bin frequencies, converted to Hz
frequency_bins = np.fft.fftfreq(number_of_samples, 1 / SAMPLING_FREQUENCY)


# ---------------- 2. Analog Butterworth prototype ----------------
# The bilinear transform warps the frequency axis, so the analog
# cutoff is prewarped:   Wc = 2*fs*tan(pi*fc/fs)
analog_cutoff = 2 * SAMPLING_FREQUENCY * np.tan(
    np.pi * CUTOFF_FREQUENCY / SAMPLING_FREQUENCY)

# The N poles sit on a circle of radius Wc, in the left half-plane
pole_index = np.arange(filter_order)
analog_poles = analog_cutoff * np.exp(
    1j * np.pi * (2 * pole_index + filter_order + 1) / (2 * filter_order))

# Gain chosen so that Ha(0) = 1  ->  no attenuation at DC
analog_gain = np.prod(-analog_poles)          # Ha(0) = 1


def filter_response(frequencies):
    """H(f) of the digital filter, obtained by the bilinear transform."""
    # On the unit circle the delay operator is z^-1 = exp(-j*w)
    z_inverse = np.exp(-2j * np.pi * frequencies / SAMPLING_FREQUENCY)

    # Bilinear transform:  s = 2*fs*(1 - z^-1)/(1 + z^-1)
    # At the Nyquist frequency z^-1 = -1, so s -> infinity and H = 0
    response = np.zeros(len(frequencies), dtype=complex)
    finite = np.abs(1 + z_inverse) > 1e-12    # z=-1 -> s=inf -> H=0
    s = 2 * SAMPLING_FREQUENCY * (1 - z_inverse[finite]) / (1 + z_inverse[finite])

    # Ha(s) = K / product(s - pole)
    response[finite] = analog_gain / np.prod(s[:, None] - analog_poles, axis=1)
    return response


# ---------------- 3. Filtering in the frequency domain ----------------
input_spectrum = np.fft.fft(input_signal)              # x[n] -> X[k]
sampled_response = filter_response(frequency_bins)     # H[k]
output_spectrum = sampled_response * input_spectrum    # multiplication
filtered_signal = np.fft.ifft(output_spectrum).real    # y[n]

# Note: multiplying two N-point DFTs is a CIRCULAR convolution,
# so the first samples of y[n] carry a wrap-around transient.


# ---------------- 4. Plots ----------------
# fftshift only reorders the bins for display, from -fs/2 to +fs/2
shifted_frequency = np.fft.fftshift(frequency_bins)
shifted_input = np.fft.fftshift(input_spectrum)
shifted_output = np.fft.fftshift(output_spectrum)

# A dense grid gives a smooth curve, independent of the signal length
dense_frequency = np.linspace(1, SAMPLING_FREQUENCY / 2, 1000)
response_db = 20 * np.log10(np.maximum(np.abs(filter_response(dense_frequency)), 1e-10))

figure, axes = plt.subplots(2, 2, figsize=(13, 8))

# Top left: filter magnitude response in dB
axes[0, 0].semilogx(dense_frequency, response_db)
axes[0, 0].axvline(CUTOFF_FREQUENCY, linestyle="--", color="red")
axes[0, 0].set(title=f"Butterworth response, N={filter_order}",
               xlabel="Frequency [Hz]", ylabel="Magnitude [dB]")

# Top right: the input signal in the time domain
axes[0, 1].stem(sample_indices, input_signal, basefmt=" ")
axes[0, 1].set(title="Input x[n]", xlabel="n", ylabel="x[n]")

# Bottom left: spectrum before and after filtering
axes[1, 0].stem(shifted_frequency, np.abs(shifted_input),
                linefmt="C0-", markerfmt="C0o", basefmt=" ", label="|X[k]|")
axes[1, 0].stem(shifted_frequency, np.abs(shifted_output),
                linefmt="C1-", markerfmt="C1s", basefmt=" ", label="|Y[k]|")
axes[1, 0].set(title="Spectrum before and after filtering",
               xlabel="Frequency [Hz]", ylabel="Magnitude")
axes[1, 0].legend()

# Bottom right: the filtered signal in the time domain
axes[1, 1].stem(sample_indices, filtered_signal, basefmt=" ")
axes[1, 1].set(title="Filtered y[n]", xlabel="n", ylabel="y[n]")

for graph in axes.flat:
    graph.grid(True)

plt.tight_layout()
plt.show()
