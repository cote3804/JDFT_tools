# converting units and types

def bias_float_to_str(bias:float) -> str:
    return f"{bias:.2f}" + "V"

def bias_str_to_float(bias:str) -> float:
    return float(bias.strip('V'))

def mu_to_she(mu:float, solvent="CANDLE") -> float:
    if solvent == "CANDLE":
        return -H_to_eV(mu) - 4.66

def H_to_eV(H:float) -> float:
    return H * 27.2114