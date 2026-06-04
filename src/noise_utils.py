import _cultiv_noise as cultiv_noise
import _honey_comb_noise as hc_noise


FLEXIBLE_PAIRWISE_PROFILES = {
    "measurement-limited": {"p_m_bulk": 1.0, "p_spam": 0.000001, "p_depolarizing": 0.000001},
    "depolarizing-limited": {"p_m_bulk": 0.000001, "p_spam": 0.000001, "p_depolarizing": 1.0},
    "spam-limited": {"p_m_bulk": 0.000001, "p_spam": 1.0, "p_depolarizing": 0.000001},
    "balanced": {"p_m_bulk": 1.0, "p_spam": 1.0, "p_depolarizing": 1.0},
}


def get_noise_model(name: str, p: float, *, profile: str | None = None):
    """
    Returns a noise model based on the given name.
    For 'flexible_pairwise', use profile to select a named ratio preset.
    """
    match name:
        case "em3" | "EM3":
            return hc_noise.NoiseModel.EM3_v3(p)
        case "SIEM3000" | "siem3000":
            return hc_noise.NoiseModel.SIEM3000(p)
        case "uniform":
            return cultiv_noise.NoiseModel.uniform_depolarizing(p)
        case "normalized_pairwise_errors":
            return cultiv_noise.NoiseModel.normalized_pairwise_errors(p)
        case "phenomenological":
            return cultiv_noise.NoiseModel.uniform_depolarizing(p, single_qubit_only=True)
        case "si1000":
            return cultiv_noise.NoiseModel.si1000(p)
        case "flexible_pairwise":
            if profile is None:
                profile = "balanced"
            ratios = FLEXIBLE_PAIRWISE_PROFILES[profile]
            return cultiv_noise.NoiseModel.flexible_pairwise(
                p_m_bulk=p * ratios["p_m_bulk"],
                p_spam=p * ratios["p_spam"],
                p_depolarizing=p * ratios["p_depolarizing"],
            )
        case _:
            raise ValueError(f"Unknown noise model name: {name}")
