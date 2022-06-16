"""Get resonance data from ImagingReso library."""

import csv
import pickle
import ImagingReso
from ImagingReso.resonance import Resonance
from scipy.signal import find_peaks


def read_elements(
    folder=ImagingReso.__path__[0],
    fname="/reference_data/ENDF_VII/_elements_list.csv",
):
    """Get available elements in ImagingReso."""
    with open(folder + fname) as f:
        reader = csv.reader(f)
        elements = [x[1] for x in reader][1:]
    return elements


def get_sigmas(elements, E_min=0.1, E_max=200.0, E_res=0.01, d=1):
    """Get abundance and cross section data for specified elements."""
    o_reso = Resonance(energy_min=E_min, energy_max=E_max, energy_step=E_res)
    # add all elements to Resonance stack
    for i in elements:
        _layer = i
        _thickness = d
        print(f"Getting resonances for {i}")
        try:
            o_reso.add_layer(formula=_layer, thickness=_thickness)
        except TypeError:
            print(
                f"Issue retrieving cross section for element: {i}",
                "\n(likely issue: no known density)",
            )
        except ValueError:
            print(
                f"Issue retrieving cross section for element: {i}",
                "\n(likely issue: no known natural isotopic ratio)",
            )
    return o_reso.stack_sigma


def save_pkl(data, fname="sigmas.pkl"):
    """Save data to pickle."""
    with open(fname, "wb") as f:
        pickle.dump([data], f)


def get_resos(sigmas, height=50, prominence=25):
    """Get all resonances in energy range above threshold."""
    dict_resos = {}
    for key in sigmas:
        isos = [
            iso for iso in sorted(sigmas[key][key].keys()) if iso[0].isdigit()
        ]
        for i, iso in enumerate(isos):
            xs = sigmas[key][key][iso]["sigma_b_raw"]
            peaks, properties = find_peaks(
                xs, height=height, prominence=prominence
            )
            peak_energies = [
                sigmas[key][key][iso]["energy_eV"][i] for i in peaks
            ]
            for j, peak in enumerate(peaks):
                if j == 0:
                    dict_resos[iso] = {}
                e_str = str(peak_energies[j])
                dict_resos[iso][e_str] = {
                    "Abundance": sigmas[key][key]["isotopic_ratio"][i],
                    "Energy": peak_energies[j],
                    "Sigma": properties["peak_heights"][j],
                }
    return dict_resos


if __name__ == "__main__":
    # Generate resonance parameters for all elements
    if "Sigma" not in locals():
        # check for pickled data
        try:
            with open("sigmas.pkl", "rb") as f:
                sigmas = pickle.load(f)
        except FileNotFoundError:
            print(
                "Did not find pickled cross section data.",
                "\nReading cross section data from ImagingReso library...",
            )
            elements = read_elements()
            sigmas = get_sigmas(elements)

    # find elements with cross-sections above specified threshold
    resos = get_resos(sigmas)

    # optional: pickle the cross section and resonance data
    pkl_flag = False
    if pkl_flag:
        save_pkl(sigmas, fname="sigmas.pkl")
        save_pkl(resos, fname="resos.pkl")
