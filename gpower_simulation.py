import sys
import json
import math
from scipy import stats

VERSION_INFO = {
    "contract_version": "1.0.0",
    "name": "Kozmos-GPowerEngine",
    "engine": "Python 3.10 + SciPy"
}

def error_exit(msg, code="INVALID_VALUE", field=None, exit_code=2):
    err_dict = {"error": msg, "code": code}
    if field:
        err_dict["field"] = field
    print(json.dumps(err_dict))
    sys.exit(exit_code)

def parse_args():
    args = sys.argv[1:]
    if "--version" in args:
        print(json.dumps(VERSION_INFO))
        sys.exit(0)
        
    parsed = {}
    i = 0
    while i < len(args):
        arg = args[i]
        if arg.startswith("--"):
            key = arg[2:]
            if i + 1 < len(args) and not args[i+1].startswith("--"):
                val = args[i+1]
                try:
                    if "." in val:
                        parsed[key] = float(val)
                    else:
                        parsed[key] = int(val)
                except ValueError:
                    parsed[key] = val
                i += 2
            else:
                parsed[key] = True
                i += 1
        else:
            i += 1
    return parsed

def validate_ranges(params):
    if "tails" in params and params["tails"] not in [1, 2]:
        error_exit("Tails 1 veya 2 olmalıdır", "OUT_OF_RANGE", "tails")
    if "alpha" in params and not (0 < params["alpha"] < 1):
        error_exit("Alpha 0 ile 1 arasında olmalıdır", "OUT_OF_RANGE", "alpha")
    if "power" in params and not (0 < params["power"] < 1):
        error_exit("Power 0 ile 1 arasında olmalıdır", "OUT_OF_RANGE", "power")
    for key in ["d", "dz", "f", "w", "f2", "ratio"]:
        if key in params and params[key] <= 0:
            error_exit(f"{key} 0'dan büyük olmalıdır", "OUT_OF_RANGE", key)
    for key in ["n1", "n2", "n_total", "groups", "df", "predictors"]:
        if key in params and params[key] < 2 and key != "df":
            error_exit(f"{key} 2 veya daha büyük olmalıdır", "OUT_OF_RANGE", key)

def calculate_t2_indep(params, mode):
    if mode == "a-priori":
        for r in ["tails", "d", "alpha", "power", "ratio"]:
            if r not in params:
                error_exit(f"Eksik parametre: {r}", "MISSING_PARAM", r)
        tails, d, alpha, target_power, ratio = int(params["tails"]), float(params["d"]), float(params["alpha"]), float(params["power"]), float(params["ratio"])
        
        n1 = 2
        while n1 < 100000:
            n2 = int(round(n1 * ratio))
            df = n1 + n2 - 2
            t_crit = stats.t.ppf(1 - alpha / tails, df)
            ncp = d * math.sqrt((n1 * n2) / (n1 + n2))
            power = 1 - (stats.nct.cdf(t_crit, df, ncp) - stats.nct.cdf(-t_crit, df, ncp)) if tails == 2 else 1 - stats.nct.cdf(t_crit, df, ncp)
            
            if power >= target_power:
                return {
                    "ncp": round(ncp, 6),
                    "critical_t": round(t_crit, 6),
                    "df": int(df),
                    "n1": int(n1),
                    "n2": int(n2),
                    "n_total": int(n1 + n2),
                    "actual_power": round(power, 6)
                }
            n1 += 1
        error_exit("Çözüm bulunamadı", "NO_SOLUTION", exit_code=1)
    elif mode == "post-hoc":
        for r in ["tails", "d", "alpha", "n1", "n2"]:
            if r not in params:
                error_exit(f"Eksik parametre: {r}", "MISSING_PARAM", r)
        tails, d, alpha, n1, n2 = int(params["tails"]), float(params["d"]), float(params["alpha"]), int(params["n1"]), int(params["n2"])
        df = n1 + n2 - 2
        t_crit = stats.t.ppf(1 - alpha / tails, df)
        ncp = d * math.sqrt((n1 * n2) / (n1 + n2))
        power = 1 - (stats.nct.cdf(t_crit, df, ncp) - stats.nct.cdf(-t_crit, df, ncp)) if tails == 2 else 1 - stats.nct.cdf(t_crit, df, ncp)
        return {
            "ncp": round(ncp, 6),
            "critical_t": round(t_crit, 6),
            "df": int(df),
            "power": round(power, 6)
        }
    else:
        error_exit(f"Desteklenmeyen mod: {mode}", "UNSUPPORTED")

def calculate_t_paired(params, mode):
    if mode == "a-priori":
        for r in ["tails", "dz", "alpha", "power"]:
            if r not in params:
                error_exit(f"Eksik parametre: {r}", "MISSING_PARAM", r)
        tails, dz, alpha, target_power = int(params["tails"]), float(params["dz"]), float(params["alpha"]), float(params["power"])
        
        n = 2
        while n < 100000:
            df = n - 1
            t_crit = stats.t.ppf(1 - alpha / tails, df)
            ncp = dz * math.sqrt(n)
            power = 1 - (stats.nct.cdf(t_crit, df, ncp) - stats.nct.cdf(-t_crit, df, ncp)) if tails == 2 else 1 - stats.nct.cdf(t_crit, df, ncp)
            
            if power >= target_power:
                return {
                    "ncp": round(ncp, 6),
                    "critical_t": round(t_crit, 6),
                    "df": int(df),
                    "n_total": int(n),
                    "actual_power": round(power, 6)
                }
            n += 1
        error_exit("Çözüm bulunamadı", "NO_SOLUTION", exit_code=1)
    elif mode == "post-hoc":
        for r in ["tails", "dz", "alpha", "n_total"]:
            if r not in params:
                error_exit(f"Eksik parametre: {r}", "MISSING_PARAM", r)
        tails, dz, alpha, n = int(params["tails"]), float(params["dz"]), float(params["alpha"]), int(params["n_total"])
        df = n - 1
        t_crit = stats.t.ppf(1 - alpha / tails, df)
        ncp = dz * math.sqrt(n)
        power = 1 - (stats.nct.cdf(t_crit, df, ncp) - stats.nct.cdf(-t_crit, df, ncp)) if tails == 2 else 1 - stats.nct.cdf(t_crit, df, ncp)
        return {
            "ncp": round(ncp, 6),
            "critical_t": round(t_crit, 6),
            "df": int(df),
            "power": round(power, 6)
        }
    else:
        error_exit(f"Desteklenmeyen mod: {mode}", "UNSUPPORTED")

def calculate_anova_oneway(params, mode):
    if mode == "a-priori":
        for r in ["f", "alpha", "power", "groups"]:
            if r not in params:
                error_exit(f"Eksik parametre: {r}", "MISSING_PARAM", r)
        f_effect, alpha, target_power, groups = float(params["f"]), float(params["alpha"]), float(params["power"]), int(params["groups"])
        
        n = groups * 2
        while n < 200000:
            df1 = groups - 1
            df2 = n - groups
            if df2 <= 0:
                n += 1
                continue
            f_crit = stats.f.ppf(1 - alpha, df1, df2)
            ncp = (f_effect ** 2) * n
            power = 1 - stats.ncf.cdf(f_crit, df1, df2, ncp)
            
            if power >= target_power:
                return {
                    "ncp": round(ncp, 6),
                    "critical_f": round(f_crit, 6),
                    "df1": int(df1),
                    "df2": int(df2),
                    "n_total": int(n),
                    "actual_power": round(power, 6)
                }
            n += 1
        error_exit("Çözüm bulunamadı", "NO_SOLUTION", exit_code=1)
    elif mode == "post-hoc":
        for r in ["f", "alpha", "n_total", "groups"]:
            if r not in params:
                error_exit(f"Eksik parametre: {r}", "MISSING_PARAM", r)
        f_effect, alpha, n, groups = float(params["f"]), float(params["alpha"]), int(params["n_total"]), int(params["groups"])
        df1 = groups - 1
        df2 = n - groups
        f_crit = stats.f.ppf(1 - alpha, df1, df2)
        ncp = (f_effect ** 2) * n
        power = 1 - stats.ncf.cdf(f_crit, df1, df2, ncp)
        return {
            "ncp": round(ncp, 6),
            "critical_f": round(f_crit, 6),
            "df1": int(df1),
            "df2": int(df2),
            "power": round(power, 6)
        }
    else:
        error_exit(f"Desteklenmeyen mod: {mode}", "UNSUPPORTED")

def calculate_chi2_contingency(params, mode):
    if mode == "a-priori":
        for r in ["w", "alpha", "power", "df"]:
            if r not in params:
                error_exit(f"Eksik parametre: {r}", "MISSING_PARAM", r)
        w, alpha, target_power, df = float(params["w"]), float(params["alpha"]), float(params["power"]), int(params["df"])
        
        n = 5
        while n < 200000:
            chi2_crit = stats.chi2.ppf(1 - alpha, df)
            ncp = (w ** 2) * n
            power = 1 - stats.ncx2.cdf(chi2_crit, df, ncp)
            
            if power >= target_power:
                return {
                    "ncp": round(ncp, 6),
                    "critical_chi2": round(chi2_crit, 6),
                    "df": int(df),
                    "n_total": int(n),
                    "actual_power": round(power, 6)
                }
            n += 1
        error_exit("Çözüm bulunamadı", "NO_SOLUTION", exit_code=1)
    elif mode == "post-hoc":
        for r in ["w", "alpha", "n_total", "df"]:
            if r not in params:
                error_exit(f"Eksik parametre: {r}", "MISSING_PARAM", r)
        w, alpha, n, df = float(params["w"]), float(params["alpha"]), int(params["n_total"]), int(params["df"])
        chi2_crit = stats.chi2.ppf(1 - alpha, df)
        ncp = (w ** 2) * n
        power = 1 - stats.ncx2.cdf(chi2_crit, df, ncp)
        return {
            "ncp": round(ncp, 6),
            "critical_chi2": round(chi2_crit, 6),
            "df": int(df),
            "power": round(power, 6)
        }
    else:
        error_exit(f"Desteklenmeyen mod: {mode}", "UNSUPPORTED")

def calculate_multireg_r2(params, mode):
    if mode == "a-priori":
        for r in ["f2", "alpha", "power", "predictors"]:
            if r not in params:
                error_exit(f"Eksik parametre: {r}", "MISSING_PARAM", r)
        f2, alpha, target_power, predictors = float(params["f2"]), float(params["alpha"]), float(params["power"]), int(params["predictors"])
        
        n = predictors + 5
        while n < 200000:
            df1 = predictors
            df2 = n - predictors - 1
            if df2 <= 0:
                n += 1
                continue
            f_crit = stats.f.ppf(1 - alpha, df1, df2)
            ncp = f2 * n
            power = 1 - stats.ncf.cdf(f_crit, df1, df2, ncp)
            
            if power >= target_power:
                return {
                    "ncp": round(ncp, 6),
                    "critical_f": round(f_crit, 6),
                    "df1": int(df1),
                    "df2": int(df2),
                    "n_total": int(n),
                    "actual_power": round(power, 6)
                }
            n += 1
        error_exit("Çözüm bulunamadı", "NO_SOLUTION", exit_code=1)
    elif mode == "post-hoc":
        for r in ["f2", "alpha", "n_total", "predictors"]:
            if r not in params:
                error_exit(f"Eksik parametre: {r}", "MISSING_PARAM", r)
        f2, alpha, n, predictors = float(params["f2"]), float(params["alpha"]), int(params["n_total"]), int(params["predictors"])
        df1 = predictors
        df2 = n - predictors - 1
        f_crit = stats.f.ppf(1 - alpha, df1, df2)
        ncp = f2 * n
        power = 1 - stats.ncf.cdf(f_crit, df1, df2, ncp)
        return {
            "ncp": round(ncp, 6),
            "critical_f": round(f_crit, 6),
            "df1": int(df1),
            "df2": int(df2),
            "power": round(power, 6)
        }
    else:
        error_exit(f"Desteklenmeyen mod: {mode}", "UNSUPPORTED")

def main():
    params = parse_args()
    test = params.get("test")
    mode = params.get("mode")
    
    if not test:
        error_exit("Test kimliği belirtilmedi (--test)", "MISSING_PARAM", "test")
    if not mode:
        error_exit("Analiz modu belirtilmedi (--mode)", "MISSING_PARAM", "mode")
        
    validate_ranges(params)
    
    if test == "t2-indep":
        result = calculate_t2_indep(params, mode)
    elif test == "t-paired":
        result = calculate_t_paired(params, mode)
    elif test == "anova-oneway":
        result = calculate_anova_oneway(params, mode)
    elif test == "chi2-contingency":
        result = calculate_chi2_contingency(params, mode)
    elif test == "multireg-r2":
        result = calculate_multireg_r2(params, mode)
    else:
        error_exit(f"Bilinmeyen test: {test}", "UNSUPPORTED", "test")
        
    print(json.dumps(result))
    sys.exit(0)

if __name__ == "__main__":
    main()