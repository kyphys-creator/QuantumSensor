from pathlib import Path
import numpy as np

DATA_DIR = Path(__file__).resolve().parents[2] / "data"


def load_eta(eta: str, mass: str) -> np.ndarray:
    return np.genfromtxt(DATA_DIR / "Eta_data" / f"Eta{eta}M{mass}_Ko.csv", delimiter=",")


def load_response_matrix(material: str, q: str, mass: str) -> np.ndarray:
    path = DATA_DIR / "curlyRplotting" / f"{material}_R10_q{q}M{mass}_Ko_thresholded_refined_exposure.csv"
    return np.genfromtxt(path, delimiter=",")


def load_ratebin(eta: str, material: str, q: str, mass: str) -> np.ndarray:
    path = DATA_DIR / "Ratebin" / f"Event10Eta{eta}Mat{material}q{q}M{mass}.csv"
    return np.genfromtxt(path, delimiter=",")
